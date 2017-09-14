/* Main file and entry point.
 *	C implementation of a neural network simulator developed during the CEMRACS 2017.
 *	This file contains the main function of the simulator.
 *
 *	Neuons are modeled as PDE with partial equation: V_t^i = V_0^i + \int_0^t b(V_S^i)dS + \sigma*W_t^i + \sum_{j=1}^N J^{j->i}*M_t^j - M_t^i * (S^i-V^{1,i}).
 *	For simplification, b, which must be lipschitzian, is taken equal to x->b(x)=-\lambda * (x - a), with '\lambda' and 'a' that depend on the neuron.
 */
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <limits.h>

#include "./SFMT-src-1.5.1/SFMT.h"		/* Simple precision rng */
#include "./dSFMT-src-2.2.3/dSFMT.h"	/* Double precision rng */

/* DEFINES */
#define 	MAX_MEMORY	100000000
#define 	BOOL		unsigned char
#define 	EPSILON		10e-7

typedef struct {
	long double lower_bound;
	long double upper_bound;
} time_interval_t;

/* GLOBAL VARIABLES */
unsigned int	spiking_times_array_index,
				nb_neurons,	/* Denoted N hereafter */
				nb_threads,
				nb_itr;

long double	conn_prob,	/* Probability of a connection between two neurons, possibly self. */
			max_time,
/* Probability of spike function */
			slope,		/*  */
			exponent,	/*  */
			time_step,	/*  */
/* Max variance for Brownian motion */
			sigma,
			sigma_squared,	/*  */
			mult_const;		/*  */

BOOL	use_rec;

/* The rngs */
sfmt_t	*	int_rngs;		/* Array [1..nb_threads+1]: one per thread plus one for the main thread */
dsfmt_t	*	double_rngs;	/* Array [1..nb_threads+1]: one per thread plus one for the main thread */

/* Result files */
FILE* f_results;

/* The b function */
long double 	*	lambda,	/* Array [1..N]: lambda constant values */
				*	a,		/* Array [1..N]: constant value at origin for the b function (see header). */
/* The spiking times */
				*	spiking_times,	/* Array [1..N]: the times of true spikes. The array is written down the result file when full. */
/* Potentials */
				*	reset_value,	/* Array [1..N]: reset values after spiking */
				*	indices,		/* Array [1..N]: numbers from 0 to nb_neurons-1, used for initialising the coupling graph. */
				*	t_last,			/* Array [1..N]: last spiking times */
				*	t_last_true,	/* Array [1..N]: last true spiking times */
				*	y_last,			/* Array [1..N]: potentials at last spiking time */
				*	var,			/* Array [1..N]: precomputed constant part for the brownian motion; sigma^2 / (2 * lambda_{i}), i the neuron index */
/* Maxima of probability of spike function */
				*	threshold,		/* Array [1..N]: threshold in the probability function */
				*	max,			/* Array [1..N]: maximum of probability function between [T_{i-1}, T_i] */
				*	max_cum_sum;	/* Array [1..N]: cumulative sum of probability function between [T_{i-1}, T_i] */
/* Interactions */
unsigned int	*	nb_couplings;			/* Array [1..N]: number of postsynaptic neurons */
sfmt_t			*	reconstruction_graph;	/* Array [1..N]: the state of the rng for reconstruction */
unsigned int	**	interaction_graph;		/* Matrix [1..N][0..nb_couplings]: for each neuron, the indices of the post-synaptic neurons or the indices of unconnected neurons, depending on the total number of neurons */

/* FUNCTIONS PROTOTYPES */
void check_alloc( const void* const var, const char * const var_name );
void init(void);
void end(void);
void compute_m_i_s( const time_interval_t* const time_intn );
void simulate(void);
long double interaction(void);
long double probability( const unsigned int neuron );
long double max_prob( const unsigned int neuron, const time_interval_t* const time_int );
unsigned int which_spiking(void);
unsigned int get_int( const unsigned int limit );
BOOL is_real_spike( const unsigned int spiking_neuron );

/* FUNCTIONS IMPLEMENTATION */
int main( int argc, char** const argv ) {
	if( argc == 4 ) {
		nb_neurons = strtoul( argv[ 1 ], NULL, 10 );
		conn_prob = strtod( argv[ 2 ], NULL );
		nb_threads = strtoul( argv[ 3 ], NULL, 10 );
	} else {
		nb_neurons = 1;
		conn_prob = 0.0;
		nb_threads = 1;
	}
	max_time	= 10;
	nb_itr		= 1.5 * nb_neurons;

	init();

	fprintf( stdout, "Hello World! There are %i arguments.\n", argc );
	fprintf( stdout, "Number of neurons: %u.\n", nb_neurons );
	fprintf( stdout, "Probability of a connection: %Lf.\n", conn_prob );

	f_results = fopen( "result.bin", "wb" );
	simulate();
	if( spiking_times_array_index > 0 )
		fwrite( spiking_times, sizeof( long double ), spiking_times_array_index, f_results );
	fclose( f_results );

	end();

	return 0;
}

inline
void check_alloc( const void* const var, const char * const var_name ) {
	if( var == NULL ) fprintf( stderr, "Cannot allocate variable %s\n", var_name );
}

void init(void) {
	unsigned int i, j, max_ind, ind;

	spiking_times_array_index = 0;
	use_rec = conn_prob * nb_neurons * nb_neurons > MAX_MEMORY && (1 - conn_prob) * nb_neurons * nb_neurons > MAX_MEMORY;
	time_step = 1e-3;
	sigma = 1.0;
	sigma_squared = sigma * sigma;
	mult_const = 5;
	slope = 1e5;
	exponent = 5;

	fflush(stdout);
	check_alloc( int_rngs = malloc( sizeof( sfmt_t ) * (nb_threads + 1) ), "integer rngs" );
	for( i = 0; i <= nb_threads; ++i ) {
		/* TO-DO: initialize all generators with different seeds that guarrantee non-overlapping */
		sfmt_init_gen_rand( &(int_rngs[ i ]), 12345 + i );
	}
	check_alloc( double_rngs = malloc( sizeof( dsfmt_t ) * (nb_threads + 1) ), "double rngs" );
	for( i = 0; i <= nb_threads; ++i ) {
		/* TO-DO: initialize all generators with different seeds that guarrantee non-overlapping */
		dsfmt_init_gen_rand( &(double_rngs[ i ]), 12345 + i );
	}
	check_alloc( lambda = malloc( sizeof( long double ) * nb_neurons ), "lambda" );
	for( i = 0; i < nb_neurons; ++i ) {
		lambda[ i ] = 10;
	}
	check_alloc( a = malloc( sizeof( long double ) * nb_neurons ), "a" );
	for( i = 0; i < nb_neurons; ++i ) {
		a[ i ] = 1.1;
	}
	check_alloc( spiking_times = malloc( sizeof( long double ) * nb_neurons ), "spiking_times" );
	for( i = 0; i < nb_neurons; ++i ) {
		spiking_times[ i ] = -1.0;
	}
	check_alloc( reset_value = malloc( sizeof( long double ) * nb_neurons ), "reset_value" );
	for( i = 0; i < nb_neurons; ++i ) {
		reset_value[ i ] = 0.0;
	}
	check_alloc( t_last = malloc( sizeof( long double ) * nb_neurons ), "t_last" );
	for( i = 0; i < nb_neurons; ++i ) {
		t_last[ i ] = 0.0;
	}
	check_alloc( t_last_true = malloc( sizeof( long double ) * nb_neurons ), "t_last_true" );
	for( i = 0; i < nb_neurons; ++i ) {
		t_last_true[ i ] = 0.0;
	}
	check_alloc( y_last = malloc( sizeof( long double ) * nb_neurons ), "y_last" );
	for( i = 0; i < nb_neurons; ++i ) {
		y_last[ i ] = 0.0;
	}
	check_alloc( var = malloc( sizeof( long double ) * nb_neurons ), "var" );
	for( i = 0; i < nb_neurons; ++i ) {
		var[ i ] = mult_const * sqrt( sigma_squared / (2 * lambda[ i ]) );
	}

	check_alloc( threshold = malloc( sizeof( long double ) * nb_neurons ), "threshold" );
	for( i = 0; i < nb_neurons; ++i ) {
		threshold[ i ] = 1.0;
	}
	check_alloc( max = malloc( sizeof( long double ) * nb_neurons ), "max" );
	for( i = 0; i < nb_neurons; ++i ) {
		max[ i ] = 0.0;
	}
	check_alloc( max_cum_sum = malloc( sizeof( long double ) * nb_neurons ), "max_cum_sum" );
	for( i = 0; i < nb_neurons; ++i ) {
		max_cum_sum[ i ] = 0.0;
	}

	check_alloc( indices = malloc( sizeof( long double ) * nb_neurons ), "indices" );
	for( i = 0; i < nb_neurons; ++i ) {
		indices[ i ] = i;
	}
	/* Draw a random integer between 0 and nb_neurons included */
	check_alloc( nb_couplings = malloc( sizeof( unsigned int ) * nb_neurons ), "nb_couplings" );
	for( i = 0; i < nb_neurons; ++i ) {
		nb_couplings[ i ] = 0;
		for( j = 0; j < nb_neurons; ++j ) {
			if( dsfmt_genrand_close1_open2( double_rngs ) < conn_prob ) {
				++(nb_couplings[ i ]);
			}
		}
	}
	if( use_rec ) {
		interaction_graph = NULL;
		check_alloc( reconstruction_graph = malloc( sizeof( sfmt_t ) * nb_neurons ), "reconstruction_graph" );
		for( i = 0; i < nb_neurons; ++i ) {
			reconstruction_graph[ i ] = *int_rngs;
			max_ind = nb_neurons - 1;
			for ( j = 0; j < nb_couplings[ i ]; ++j ) {
				get_int( max_ind );
				--max_ind;
			}
		}
	} else {
		reconstruction_graph = NULL;
		check_alloc( interaction_graph = malloc( sizeof( unsigned int* ) * nb_neurons ), "interaction_graph" );
		for( i = 0; i < nb_neurons; ++i ) {
			check_alloc( interaction_graph[ i ] = malloc( sizeof( unsigned int ) * nb_couplings[ i ] ), "interaction_graph_i" );
			max_ind = nb_neurons - 1;
			for ( j = 0; j < nb_couplings[ i ]; ++j ) {
				ind = get_int( max_ind );
				interaction_graph[ i ][ j ] = indices[ ind ];
				indices[ ind ] = indices[ max_ind ];
				indices[ max_ind ] = interaction_graph[ i ][ j ];
				--max_ind;
			}
		}
	}
}

inline
unsigned int get_int( const unsigned int limit ) {
	unsigned int ret_val, u, divide;

	divide = UINT_MAX / (limit + 1);	/* {0, ..., limit}, so limit+1 possibilities */
	do {
		u = sfmt_genrand_uint64( int_rngs );
		ret_val = u / divide;
	} while( ret_val > limit );
	return ret_val;
}

void end(void) {
	unsigned int i;

	use_rec = conn_prob * nb_neurons * nb_neurons > MAX_MEMORY;

	free( int_rngs );
	free( double_rngs );

	free( lambda );
	free( a );

	free( reset_value );
	free( t_last );
	free( t_last_true );
	free( y_last );
	free( var );

	free( threshold );
	free( max );
	free( max_cum_sum );
	
	free( nb_couplings );

	if( use_rec ) {
		free(reconstruction_graph);
	} else {
		for( i = 0; i < nb_neurons; ++i ) {
			free( interaction_graph[ i ] );
		}
		free( interaction_graph );
	}
}

inline
void compute_m_i_s( const time_interval_t* const time_int ) {
	/* Declarations */
	unsigned int i;

	/* Initializations */
	i = 0;

	/* Algorithm */
	max[ 0 ] = max_prob( 0, time_int );
	max_cum_sum[ 0 ] = max[ 0 ];
	for( i = 1; i < nb_neurons; ++i ) {
		max[ i ] = max_prob( i, time_int );
		max_cum_sum[ i ] = max_cum_sum[ i - 1 ] + max[ i ];
	}
}

void simulate(void) {
	/* Declarations */
	sfmt_t			save_rng;
	long double		delta_t,
					t,
					v;
	long double		u[2],
					normals[2];
	time_interval_t	time_int;
	unsigned int	i, ind_buff,
					max_ind,
					spiking_neuron,	/* The index of the current spiking neuron */
					itr,			/* The number of iterations (number of true spikes) since the beginning */
					normals_ind;	/* The index of the normal number used for the computation of the brownian motion. @see #normals */

	/* Initializations */
	u[ 0 ] = dsfmt_genrand_open_open( double_rngs ); u[ 1 ] = dsfmt_genrand_open_open( double_rngs );
	normals[ 0 ] = sqrt( -2 * log( u[ 0 ] ) ) * cos( u[ 1 ] ); normals[ 1 ] = sqrt( -2 * log( u[ 0 ] ) ) * sin( u[ 1 ] );

	time_int.lower_bound = -time_step;
	time_int.upper_bound = 0.0;
	
	itr = 0;
	normals_ind = 0;

	/* Algorithm */
NEXT_STEP:
/**	fprintf( f_results, "%Lf: NEXT_STEP\n", time_int.upper_bound );	**/
	delta_t = 0.0;
	do {
		time_int.lower_bound = time_int.upper_bound;
		time_int.upper_bound += (time_int.lower_bound + time_step > max_time)? (max_time - time_int.lower_bound): time_step;
		if( time_int.upper_bound >= max_time || itr > nb_itr ) return;
		compute_m_i_s( &time_int );
	} while( fabs( max_cum_sum[ nb_neurons - 1 ] ) < EPSILON );
	/* TO-DO: parallellize the following */
NEXT_SPIKE:
	v = -log( dsfmt_genrand_open_close( double_rngs ) );
	delta_t += v / max_cum_sum[ nb_neurons - 1 ];	/* Generating a time step of exponential distribution E(\sum max_i) */
	if( time_int.lower_bound + delta_t > time_int.upper_bound || itr > nb_itr ) goto NEXT_STEP;	/* If next potential spike is over time bound, must change of time interval */
/**	fprintf( f_results, "In bound\n" );	**/
	spiking_neuron = which_spiking();
	/* TO-DO: update the potential of spiking neuron before computing if the spike is real. */
	t = time_int.lower_bound + delta_t - t_last[ spiking_neuron ];
	t_last[ spiking_neuron ]	= time_int.lower_bound + delta_t;
	y_last[ spiking_neuron ]	= a[ spiking_neuron ]
									+ exp( -lambda[ spiking_neuron ] * t ) * (y_last[ spiking_neuron ] - a[ spiking_neuron ])
									+ sqrt( sigma_squared * (1 - exp( -2 * lambda[ spiking_neuron ] * t )) / (2 * lambda[ spiking_neuron ]) ) * normals[ normals_ind ];

	if( normals_ind ) {
		u[ 0 ] = dsfmt_genrand_open_open( double_rngs ); u[ 1 ] = dsfmt_genrand_open_open( double_rngs );
		normals[ 0 ] = sqrt( -2 * log( u[ 0 ] ) ) * cos( u[ 1 ] ); normals[ 1 ] = sqrt( -2 * log( u[ 0 ] ) ) * sin( u[ 1 ] );
	}
	normals_ind = (normals_ind + 1) & 0x01; /* normals_ind \in {0,1}, so normals_ind := (normals_ind + 1) % 2 */
	
	if( !(probability( spiking_neuron ) > dsfmt_genrand_close_open( double_rngs ) * max[ spiking_neuron ]) ) {
		goto NEXT_SPIKE;
	}
	/* Else post-synaptic neurons' potential is updated. */
	++itr;
	/* TO-DO: update  */
	y_last[ spiking_neuron ] = reset_value[ spiking_neuron ];

	spiking_times[ spiking_times_array_index++ ] = time_int.lower_bound + delta_t;
	if( spiking_times_array_index > nb_neurons - 1 ) {
		fwrite( spiking_times, sizeof( long double ), nb_neurons, f_results );
		spiking_times_array_index = 0;
	}
	if( use_rec ) {
		save_rng = *int_rngs;
		*int_rngs = reconstruction_graph[ spiking_neuron ];
		max_ind = nb_neurons - 1;
		for ( i = 0; i < nb_couplings[ spiking_neuron ]; ++i ) {
			ind = get_int( max_ind );
			y_last[ indices[ ind ] ] += interaction();

			ind_buff = indices[ ind ];
			indices[ ind ] = indices[ max_ind ];
			indices[ max_ind ] = ind_buff;
			
			--max_ind;
		}
		for( i = 0; i < nb_neurons; ++i ) {
			indices[ i ] = i;
		}
		*int_rngs = save_rng;
	} else {
		for( i = 1; i < nb_couplings[ spiking_neuron ]; ++i ) {
			y_last[ interaction_graph[ spiking_neuron ][ i ] ] += interaction();
		}
	}

	time_int.upper_bound = time_int.lower_bound + delta_t;
	goto NEXT_STEP;
}

inline
long double interaction(void) {
	return 1 / nb_neurons;
}

inline
unsigned int which_spiking(void) {
	long double value;
	unsigned int neuron;

	value = dsfmt_genrand_close_open( double_rngs ) * max_cum_sum[ nb_neurons - 1 ];
	for( neuron = 0; neuron < nb_neurons; ++neuron ) {
		if( max_cum_sum[ neuron ] >= value ) return neuron;
	}
	return -1;	/* Shall never be reached */
}

inline
long double probability( const unsigned int neuron ) {
	long double prob;
	prob = y_last[ neuron ] - threshold[ neuron ];
	return (prob > 0)? slope * pow( prob, exponent ): 0.0;
}

inline
long double max_prob( const unsigned int neuron, const time_interval_t* const time_int ) {
	long double max_prob;
	max_prob	= a[ neuron ]
					+ exp( -lambda[ neuron ] * (time_int->upper_bound - t_last[ neuron ]) ) * (y_last[ neuron ] - a[ neuron ])
					+ var[ neuron ] * sqrt(1 - exp( -2 * lambda[neuron] * (time_int->upper_bound - t_last[neuron]) ))
				- threshold[ neuron ];
	return (max_prob > 0)? slope * pow( max_prob, exponent ): 0.0;
}