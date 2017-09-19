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
#include <time.h>
#include <sys/stat.h>

#include "sim.h"
#include "./SFMT-src-1.5.1/SFMT.h"		/* Simple precision rng */
#include "./dSFMT-src-2.2.3/dSFMT.h"	/* Double precision rng */

/* DEFINES */
#define 	MAX_MEMORY	100000000
#define 	BOOL		unsigned char
#define 	EPSILON		10e-7

/* GLOBAL VARIABLES */
unsigned int	i, j,						/* Two generic buffer indices */
				ind_buff, max_ind,			/* Indices for the array of indices (meta) */
				spiking_times_array_index,	/* The current memory box available for storing a time value */
				nb_neurons,					/* Denoted N hereafter */
				size,						/* = nb_neurons - 1 */
				nb_itr;						/* Number of iterations (jumps) passed since the begining */

long double	conn_prob,	/* Probability of a connection between two neurons, possibly self. */
			max_time,	/* The maximum allowed time. */
/* Probability of spike function: slope * (Pos(potential-threshold))^exponent */
			slope,		/* The slope as defined above */
			exponent,	/* The exponent as defined above */
			time_step,	/*  */
/* Max variance for Brownian motion */
			sigma,			/*  */
			sigma_squared,	/*  */
			mult_const;		/*  */

BOOL	use_rec;

/* The rngs */
sfmt_t	*	int_rngs;		/* Array [1..nb_threads+1]: one per thread plus one for the main thread */
dsfmt_t	*	double_rngs;	/* Array [1..nb_threads+1]: one per thread plus one for the main thread */

/* Seeds, parameters and results files */
FILE		*	f_results,
			*	f_seeds,
			*	f_input,
			*	f_output;
uint32_t	*	seed;

/* The b function */
long double 	*	lambda,	/* Array [1..N]: lambda constant values */
				*	a,		/* Array [1..N]: constant value at origin for the b function (see header). */
/* The spiking times */
				*	spiking_times,	/* Array [1..N]: the times of true spikes. The array is written down the result file when full. */
/* Potentials */
				*	reset_value,	/* Array [1..N]: reset values after spiking */
				*	t_last,			/* Array [1..N]: last spiking times */
				*	t_last_true,	/* Array [1..N]: last true spiking times */
				*	y_last,			/* Array [1..N]: potentials at last spiking time */
				*	var,			/* Array [1..N]: precomputed constant part for the brownian motion; sigma^2 / (2 * lambda_{i}), i the neuron index */
/* Probability of spike function and maxima */
				*	threshold,		/* Array [1..N]: threshold in the probability function */
				*	max,			/* Array [1..N]: maximum of probability function between [T_{i-1}, T_i] */
				*	max_cum_sum;	/* Array [1..N]: cumulative sum of probability function between [T_{i-1}, T_i] */
/* Interactions */
sfmt_t			*	reconstruction_graph;	/* Array [1..N]: the state of the rng for reconstruction */
unsigned int	*	nb_couplings,			/* Array [1..N]: number of postsynaptic neurons */
				*	indices,				/* Array [1..N]: numbers from 0 to nb_neurons-1, used for initialising the coupling graph. */
				**	interaction_graph;		/* Matrix [1..N][0..nb_couplings]: for each neuron, the indices of the post-synaptic neurons or the indices of unconnected neurons, depending on the total number of neurons */
char	*	str_f_input,
		*	str_f_output,
		*	str_f_results,
		*	str_folder;

/* FUNCTIONS IMPLEMENTATION */
int main( int argc, char** const argv ) {
	if( argc == 5 ) {
		nb_neurons	= strtoul( argv[ 1 ], NULL, 10 );
		conn_prob	= strtod( argv[ 2 ], NULL );
		max_time	= strtod( argv[ 3 ], NULL );
		nb_itr		= strtoul( argv[ 4 ], NULL, 10 );
		str_f_input	= "Results/Seeds/seeds.bin";
		f_input		= fopen( str_f_input, "rb" );
	} else if( argc == 2 ) {
		str_f_input	= argv[1];
		f_input		= fopen( str_f_input, "rb" );
		if( fread( &nb_neurons, 1, sizeof( unsigned int ), f_input ) < sizeof( unsigned int ) ) {
			fprintf( stderr, "Error on reading number of neurons from file %s.\n", argv[ 1 ] );
			exit(-1);
		}
		if( fread( &conn_prob, 1, sizeof( long double ), f_input ) < sizeof( long double ) ) {
			fprintf( stderr, "Error on reading connection probability from file %s.\n", argv[ 1 ] );
			exit(-1);
		}
		if( fread( &max_time, 1, sizeof( long double ), f_input ) < sizeof( long double ) ) {
			fprintf( stderr, "Error on reading maximum time from file %s.\n", argv[ 1 ] );
			exit(-1);
		}
		if( fread( &nb_itr, 1, sizeof( unsigned int ), f_input ) < sizeof( unsigned int ) ) {
			fprintf( stderr, "Error on reading maximum number of iterations from file %s.\n", argv[ 1 ] );
			exit(-1);
		}
	} else {
		fprintf( stdout, "Parameters must be given, either by argument or in a file (of path given by argument).\nPlease relaunch using the good arguments.\n" );
		return 0;
	}

	create();
	init();

	simulate();
	if( spiking_times_array_index > 0 ) {
		fwrite( spiking_times, sizeof( long double ), spiking_times_array_index, f_results );
	}

	destroy();

	return 0;
}

inline
void check_alloc( const void* const var, const char * const var_name ) {
	if( var == NULL ) {
		fprintf( stderr, "Cannot allocate variable %s\n", var_name );
		abort();
	}
}

void create(void) {
	use_rec = conn_prob * nb_neurons * nb_neurons > MAX_MEMORY && (1 - conn_prob) * nb_neurons * nb_neurons > MAX_MEMORY;
	
	check_alloc( str_f_output = malloc( sizeof( char ) * 100 ), "str_f_output" );
	check_alloc( str_folder = malloc( sizeof( char ) * 100 ), "str_folder" );
	check_alloc( str_f_results = malloc( sizeof( char ) * 100 ), "str_f_results" );

	check_alloc( int_rngs = malloc( sizeof( sfmt_t ) ), "integer rngs" );
	check_alloc( double_rngs = malloc( sizeof( dsfmt_t ) ), "double rngs" );
	
	check_alloc( lambda = malloc( sizeof( long double ) * nb_neurons ), "lambda" );
	check_alloc( a = malloc( sizeof( long double ) * nb_neurons ), "a" );
	
	check_alloc( spiking_times = malloc( sizeof( long double ) * nb_neurons ), "spiking_times" );
	check_alloc( reset_value = malloc( sizeof( long double ) * nb_neurons ), "reset_value" );
	check_alloc( t_last = malloc( sizeof( long double ) * nb_neurons ), "t_last" );
	check_alloc( t_last_true = malloc( sizeof( long double ) * nb_neurons ), "t_last_true" );
	check_alloc( y_last = malloc( sizeof( long double ) * nb_neurons ), "y_last" );
	check_alloc( var = malloc( sizeof( long double ) * nb_neurons ), "var" );
	check_alloc( threshold = malloc( sizeof( long double ) * nb_neurons ), "threshold" );
	check_alloc( max = malloc( sizeof( long double ) * nb_neurons ), "max" );
	check_alloc( max_cum_sum = malloc( sizeof( long double ) * nb_neurons ), "max_cum_sum" );
	check_alloc( indices = malloc( sizeof( long double ) * nb_neurons ), "indices" );
	
	check_alloc( nb_couplings = malloc( sizeof( unsigned int ) * nb_neurons ), "nb_couplings" );
	
	if( use_rec ) {
		interaction_graph = NULL;
		check_alloc( reconstruction_graph = malloc( sizeof( sfmt_t ) * nb_neurons ), "reconstruction_graph" );
	} else {
		reconstruction_graph = NULL;
		check_alloc( interaction_graph = malloc( sizeof( unsigned int* ) * nb_neurons ), "interaction_graph" );
	}
}

void init(void) {
	/* time_t t = time( NULL ); */
	sprintf( str_folder, "Results/Sim%u", (unsigned int) time( NULL ) );
	sprintf( str_f_output, "%s/seeds-%u-%1.10Lf.bin", str_folder, nb_neurons, conn_prob );
	sprintf( str_f_results, "%s/result-%u-%1.10Lf.bin", str_folder, nb_neurons, conn_prob );
	mkdir( "Results", 0700 ); mkdir( str_folder, 0700 ); mkdir( "Results/Seeds", 0700 );

	f_results		= fopen( str_f_results, "wb" );
	size			= nb_neurons - 1;
	time_step		= 1e-3;
	sigma			= 1.0;
	sigma_squared	= sigma * sigma;
	mult_const		= 5;
	slope			= 1e5;
	exponent		= 5;
	spiking_times_array_index	= 0;

	/* Seeds the RNGs */
	if( f_input == NULL ) {
		sfmt_init_gen_rand( int_rngs, 12345 );
		dsfmt_init_gen_rand( double_rngs, 12345 );
	} else {
		if( fread( int_rngs, 1, sizeof( sfmt_t ), f_input ) < sizeof( sfmt_t ) ) {
			fprintf( stderr, "Error on reading the random integer generator state from file %s.\n", str_f_input );
			exit(-1);
		}
		if( fread( double_rngs, 1, sizeof( dsfmt_t ), f_input ) < sizeof( dsfmt_t ) ) {
			fprintf( stderr, "Error on reading the random real number generator state from file %s.\n", str_f_input );
			exit(-1);
		}
		fclose( f_input );
	}
	
	/* The b function */
	for( i = 0; i < nb_neurons; ++i ) {
		lambda[ i ] = 0.1;
		a[ i ] = 1.1;
		
		spiking_times[ i ] = -1.0;
	
		reset_value[ i ] = 0.0;
	
		t_last[ i ] = 0.0;
		t_last_true[ i ] = 0.0;
		y_last[ i ] = 0.0;
		
		var[ i ] = mult_const * sqrt( sigma_squared / (2 * lambda[ i ]) );
		threshold[ i ] = 1.0;
		max[ i ] = 0.0;
		max_cum_sum[ i ] = 0.0;
		
		indices[ i ] = i;
		/* Draw a random integer between 0 and nb_neurons included */
		nb_couplings[ i ] = 0;
		for( j = 0; j < nb_neurons; ++j ) {
			if( dsfmt_genrand_close1_open2( double_rngs ) < conn_prob ) {
				++(nb_couplings[ i ]);
			}
		}
	}
	
	if( use_rec ) {
		for( i = 0; i < nb_neurons; ++i ) {
			reconstruction_graph[ i ] = *int_rngs;
			max_ind = size;
			for ( j = 0; j < nb_couplings[ i ]; ++j ) {
				get_int( max_ind );
				--max_ind;
			}
		}
	} else {
		for( i = 0; i < nb_neurons; ++i ) {
			max_ind = size;
			check_alloc( interaction_graph[ i ] = malloc( sizeof( unsigned int ) * nb_couplings[ i ] ), "interaction_graph_i" );
			for ( j = 0; j < nb_couplings[ i ]; ++j ) {
				ind_buff = get_int( max_ind );
				interaction_graph[ i ][ j ] = indices[ ind_buff ];
				indices[ ind_buff ] = indices[ max_ind ];
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

void save(void) {
	f_output = fopen( str_f_output, "wb" );
	f_seeds = fopen( "Results/Seeds/seeds.bin", "wb" );
	fwrite( &nb_neurons, 1, sizeof( unsigned int ), f_output );
	fwrite( &conn_prob, 1, sizeof( long double ), f_output );
	fwrite( &max_time, 1, sizeof( long double ), f_output );
	fwrite( &nb_itr, 1, sizeof( unsigned int ), f_output );
	fwrite( int_rngs, 1, sizeof( sfmt_t ), f_output );
	fwrite( double_rngs, 1, sizeof( dsfmt_t ), f_output );
	
	fwrite( int_rngs, 1, sizeof( sfmt_t ), f_seeds );
	fwrite( double_rngs, 1, sizeof( dsfmt_t ), f_seeds );
	
	fclose( f_output );
	fclose( f_seeds );
}

void destroy(void) {
	fclose( f_results );
	
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
	sfmt_t			save_rng;		/* An RNG state, used before the reconstruction of the interaction graph */
	long double		delta_t,		/* The time advance since the lower bound of the time interval */
					t,				/* A buffer for a time value */
					u[2],			/* Two uniformly generated values, used for the genration of normaly distributed numbers */
					normals[2];		/* Two normaly distributed independent numbers */
					sqrt_log;		/* A buffer for an intermediate value for the computation of the normal numbers */
	time_interval_t	time_int;		/* A time interval for the computation of the Poisson process */
	unsigned int	spiking_neuron,	/* The index of the current spiking neuron */
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
	/* Wait for some potentials to raise high enough */
	delta_t = 0.0;
	do {
		time_int.lower_bound = time_int.upper_bound;
		time_int.upper_bound += (time_int.lower_bound + time_step > max_time)? (max_time - time_int.lower_bound): time_step;
		if( time_int.upper_bound >= max_time || itr > nb_itr ) return;
		compute_m_i_s( &time_int );
	} while( fabs( max_cum_sum[ size ] ) < EPSILON );
NEXT_SPIKE:
	delta_t += -log( dsfmt_genrand_open_close( double_rngs ) ) / max_cum_sum[ size ];			/* Generating a time step of exponential distribution E(\sum max_i) */
	if( time_int.lower_bound + delta_t > time_int.upper_bound || itr > nb_itr ) goto NEXT_STEP;	/* If next potential spike is over time bound, must change of time interval */
	spiking_neuron = which_spiking();
	t = time_int.lower_bound + delta_t - t_last[ spiking_neuron ];
	t_last[ spiking_neuron ]	= time_int.lower_bound + delta_t;
	y_last[ spiking_neuron ]	= a[ spiking_neuron ]
									+ exp( -lambda[ spiking_neuron ] * t ) * (y_last[ spiking_neuron ] - a[ spiking_neuron ])
									+ sqrt( sigma_squared * (1 - exp( -2 * lambda[ spiking_neuron ] * t )) / (2 * lambda[ spiking_neuron ]) ) * normals[ normals_ind ];

	if( normals_ind ) {
		/* Marsaglia's method for generating normally distributed independent random numbers */
		do {
			u[ 0 ]		= 2.0 * dsfmt_genrand_open_open( double_rngs ) - 1.0;	/* cos( U([0,1]) ) ~ U([-1,1]) */
			u[ 1 ]		= 2.0 * dsfmt_genrand_open_open( double_rngs ) - 1.0;	/* sin( U([0,1]) ) ~ U([-1,1]) */
			sqrt_log	= u[ 0 ] * u[ 0 ] + u[ 1 ] * u[ 1 ];					/* Only points in the unit disk are accepted */
		} while( sqrt_log >= 1.0 );												/* The rejection rate is about 1 - \pi / 4 ~ 21%. */
		sqrt_log = sqrt( -2 * log( sqrt_log ) / sqrt_log );
		normals[ 0 ] = sqrt_log * u[ 0 ]; normals[ 1 ] = sqrt_log * u[ 1 ];
	}
	normals_ind = (normals_ind + 1) & 0x01; /* normals_ind \in {0,1}, so normals_ind := (normals_ind + 1) % 2 */
	
	if( !(probability( spiking_neuron ) > dsfmt_genrand_close_open( double_rngs ) * max[ spiking_neuron ]) ) {
		goto NEXT_SPIKE;
	}
	/* From now on the jump has been accepted,  */
	++itr;
	y_last[ spiking_neuron ] = reset_value[ spiking_neuron ];

	spiking_times[ spiking_times_array_index ] = time_int.lower_bound + delta_t;
	++spiking_times_array_index;
	if( spiking_times_array_index > size ) {
		fwrite( spiking_times, sizeof( long double ), nb_neurons, f_results );
		spiking_times_array_index = 0;
	}
	if( use_rec ) {	/* Interaction graph was too big, must be reconstructed on the fly */
		save_rng = *int_rngs;
		*int_rngs = reconstruction_graph[ spiking_neuron ];
		/* At this point it is assumed the array of indices 'indices' is sorted */
		max_ind = size;
		for( i = 0; i < nb_couplings[ spiking_neuron ]; ++i ) {
			j = get_int( max_ind ); /* A random non already chosen index in the array of indices (meta) */
			y_last[ indices[ j ] ] += interaction();

			/* The chosen index is put at the end of the array so that it will not be chosen again */
			ind_buff = indices[ j ];
			indices[ j ] = indices[ max_ind ];
			indices[ max_ind ] = ind_buff;
			/* The maximum index of the never-chosen indices is decreased as there is a new chosen index */
			--max_ind;
		}
		/* Sorting the array of indices after use */
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

	value = dsfmt_genrand_close_open( double_rngs ) * max_cum_sum[ size ];
	for( i = 0; i < nb_neurons; ++i ) {
		if( max_cum_sum[ i ] >= value ) return i;
	}
	/* Shall never be reached */
	fprintf( stderr, "Error on found spiking neuron: over max. Aborting...\n" );
	exit( -1 );
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
