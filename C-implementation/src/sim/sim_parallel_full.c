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
#include <pthread.h>
#include <omp.h>

#include "sim.h"
#include "../lib/SFMT-src-1.5.1/SFMT.h"		/* Simple precision rng */
#include "../lib/dSFMT-src-2.2.3/dSFMT.h"	/* Double precision rng */

/* DEFINES */
#define 	MAX_MEMORY	100		/* The maximum memory allowed for the storing of the interaction graph. */
#define 	BOOL		unsigned char	/* Boolean type, "equivalent" to the one defined in ++c */
#define 	EPSILON		10e-7			/* The minimum comparison value for floating point numbers.
										   Any real number of absolute value below this value is considered equal to zero.
										 */
#define		PATHSIZE	100				/* The maximum size for the paths to seeds and results files. */
#define		NB_THREADS	2				/* 8 + 1 for the main thread. */

/* GLOBAL VARIABLES */
unsigned int	i, j,						/* Two generic buffer indices */
				ind_buff, max_ind,			/* Indices for the array of indices (meta) */
				spiking_times_array_index,	/* The current memory box available for storing a time value */
				nb_neurons,					/* Denoted N hereafter */
				size,						/* = nb_neurons - 1 */
				nb_itr,						/* Maximum number of jumps (sampling acceptations) allowed for the simulation */
				nb_accepted,				/* Number of accepted jumps */
				nb_rejected,				/* Number of rejected jumps */
				nb_threads,					/* Number of threads (including the main thread) in use in the simulation, depending on the connection type (value of conn_e) */
				indices_size,				/*  */
				reconstruction_graph_size,	/*  */
				int_rngs_size,				/*  */
				threads_size,				/*  */
				spiking_neuron;				/* The index of the current spiking neuron */

long double	conn_prob,	/* Probability of a connection between two neurons, possibly self */
			max_time,	/* The maximum allowed time */
/* Probability of spike function: slope * (Pos(potential-threshold))^exponent */
			slope,		/* The slope as defined above */
			exponent,	/* The exponent as defined above */
			time_step,	/* For the thinning method, the size of the time intervals in which spikes are looked for */
/* Max variance for Brownian motion */
			sigma,			/* A control for the importance of the Brownian motion */
			sigma_squared,	/* The same as above, but squared (both are used in computations, that's why) */
			mult_const;		/* A constant, caculated from, but not only, sigma_squared (or sigma depending on the situation) */
			
conn_type	conn_e;	/* The type of connection */

/* The rngs and associated threads */
sfmt_t		*	int_rngs;		/* Array [1..NB_THREADS]: one per thread plus one for the main thread */
dsfmt_t		*	double_rngs;	/* Array [1..NB_THREADS]: one per thread plus one for the main thread */
pthread_t	*	threads;		/* Array [1..NB_THREADS]: the threads used for the parallel computation of the reconstruction */

/* Seeds, parameters and results files */
FILE		*	f_results,	/*  */
			*	f_seeds,	/*  */
			*	f_input,	/*  */
			*	f_output;	/*  */
uint32_t	*	seed;		/* The base seed used to initialize the rngs */

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
sfmt_t			**	reconstruction_graph;	/* Matrix [1..NB_THREADS][1..N]: the states of the rngs at different stages of the reconstruction */
unsigned int	*	nb_couplings,			/* Array [1..N]: number of postsynaptic neurons */
				*	thread_indices,			/* Array [1..N]: index of the thread */
				**	indices,				/* Matrix [1..NB_THREADS][1..N]: numbers from 0 to nb_neurons-1, used for initialising the coupling graph. */
				**	interaction_graph;		/* Matrix [1..N][0..nb_couplings]: for each neuron, the indices of the post-synaptic neurons or the indices of unconnected neurons, depending on the total number of neurons */
char	*	str_f_input,	/*  */
		*	str_f_output,	/*  */
		*	str_f_results,	/*  */
		*	str_folder;		/* Path to the folder of results, relatively to the folder of execution */

/* FUNCTIONS IMPLEMENTATION */
int main( int argc, char** const argv ) {
	if( argc == 5 ) {
		nb_neurons	= strtoul( argv[ 1 ], NULL, 10 );
		conn_prob	= strtod( argv[ 2 ], NULL );
		max_time	= strtod( argv[ 3 ], NULL );
		nb_itr		= strtoul( argv[ 4 ], NULL, 10 );
		str_f_input	= "../Results/Seeds/seeds.bin";
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
		fprintf( stderr, "Parameters must be given, either by argument or in a file (of path given by argument).\nPlease relaunch using the good arguments.\n" );
		return 0;
	}
fprintf( stdout, "Before create and init\n" );
fflush( stdout );
	create();
	init();
fprintf( stdout, "Before simulate\n" );
fflush( stdout );
	simulate();
	if( spiking_times_array_index > 0 ) {
		fwrite( spiking_times, sizeof( long double ), spiking_times_array_index, f_results );
	}
fprintf( stdout, "Before destroy\n" );
fflush( stdout );
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
	conn_e			=	( conn_prob * nb_neurons * nb_neurons > MAX_MEMORY )?	RECONSTRUCTION:
						( conn_prob >= 1.0 - EPSILON )?							FULL:
						( conn_prob <= EPSILON )?								INDEPENDENT:
																				RANDOM;
	int_rngs_size	=	( conn_e == RECONSTRUCTION )?	1 + NB_THREADS:
														1;
	indices_size	=	( conn_e == INDEPENDENT || conn_e == COMPLETE )?	0:
						( conn_e == RECONSTRUCTION )?						NB_THREADS:
																			1;
	threads_size	=	( conn_e == RECONSTRUCTION )?	NB_THREADS:
														0;
	reconstruction_graph_size	=	threads_size;
	
	check_alloc( str_f_output	= malloc( sizeof( char ) * PATHSIZE ), "str_f_output" );
	check_alloc( str_folder		= malloc( sizeof( char ) * PATHSIZE ), "str_folder" );
	check_alloc( str_f_results	= malloc( sizeof( char ) * PATHSIZE ), "str_f_results" );

	check_alloc( int_rngs		= malloc( sizeof( sfmt_t ) * int_rngs_size ), "integer rngs" );
	check_alloc( double_rngs	= malloc( sizeof( dsfmt_t ) ), "double rngs" );
	check_alloc( threads		= malloc( sizeof( pthread_t ) * threads_size ), "threads" );
	check_alloc( thread_indices	= malloc( sizeof( unsigned int ) * threads_size ), "thread_indices" );
	
	check_alloc( lambda	= malloc( sizeof( long double ) * nb_neurons ), "lambda" );
	check_alloc( a		= malloc( sizeof( long double ) * nb_neurons ), "a" );
	
	check_alloc( spiking_times	= malloc( sizeof( long double ) * nb_neurons ), "spiking_times" );
	check_alloc( reset_value	= malloc( sizeof( long double ) * nb_neurons ), "reset_value" );
	check_alloc( t_last			= malloc( sizeof( long double ) * nb_neurons ), "t_last" );
	check_alloc( t_last_true	= malloc( sizeof( long double ) * nb_neurons ), "t_last_true" );
	check_alloc( y_last			= malloc( sizeof( long double ) * nb_neurons ), "y_last" );
	check_alloc( var			= malloc( sizeof( long double ) * nb_neurons ), "var" );
	check_alloc( threshold		= malloc( sizeof( long double ) * nb_neurons ), "threshold" );
	check_alloc( max			= malloc( sizeof( long double ) * nb_neurons ), "max" );
	check_alloc( max_cum_sum	= malloc( sizeof( long double ) * nb_neurons ), "max_cum_sum" );
	
	check_alloc( indices = malloc( sizeof( long double* ) * indices_size ), "indices" );
	for( i = 0; i < indices_size; ++i ) {	/* One per thread */
		check_alloc( indices[ i ] = malloc( sizeof( long double ) * nb_neurons ), "indices_i" );
	}
	
	check_alloc( nb_couplings = malloc( sizeof( unsigned int ) * nb_neurons ), "nb_couplings" );
	
	switch( conn_e ) {
		case FULL:				/*  */
		case COMPLETE:			/*  */
		case INDEPENDENT:		/*  */
								interaction_graph		= NULL;
								reconstruction_graph	= NULL;
			break;
		case RECONSTRUCTION:	/*  */
								interaction_graph		= NULL;
								check_alloc( reconstruction_graph = malloc( sizeof( sfmt_t* ) * reconstruction_graph_size ), "reconstruction_graph" );
								for( i = 0; i < reconstruction_graph_size; ++i ) {
									check_alloc( reconstruction_graph[ i ] = malloc( sizeof( sfmt_t ) ), "reconstruction_graph_i" );
								}
			break;
		case RANDOM:			/*  */
								reconstruction_graph	= NULL;
								check_alloc( interaction_graph = malloc( sizeof( unsigned int* ) * nb_neurons ), "interaction_graph" );
			break;
	}
}

void init(void) {
fprintf( stdout, "In init\n" );
fflush( stdout );
	sprintf( str_folder, "../Results/Sim%u", (unsigned int) time( NULL ) );
	sprintf( str_f_output, "%s/seeds-%u-%1.10Lf.bin", str_folder, nb_neurons, conn_prob );
	sprintf( str_f_results, "%s/result-%u-%1.10Lf.bin", str_folder, nb_neurons, conn_prob );
	mkdir( "../Results", 0700 ); mkdir( str_folder, 0700 ); mkdir( "../Results/Seeds", 0700 );

	f_results		= fopen( str_f_results, "wb" );
	size			= nb_neurons - 1;
	nb_accepted		= 0;
	nb_rejected		= 0;
	time_step		= 1e-3;
	sigma			= 1.0;
	sigma_squared	= sigma * sigma;
	mult_const		= 5;
	slope			= 1e5;
	exponent		= 5;
	spiking_times_array_index	= 0;

	/* Seeds the RNGs */
	if( f_input == NULL ) {
		for( i = 0; i < nb_threads; ++i ) {
			sfmt_init_gen_rand( &(int_rngs[ i ]), 12345 + i );
		}
		dsfmt_init_gen_rand( double_rngs, 12345 );
	} else {
		if( fread( int_rngs, 1, sizeof( sfmt_t ) * int_rngs_size, f_input ) < sizeof( sfmt_t ) ) {
			fprintf( stderr, "Error on reading the random integer generator state from file %s.\n", str_f_input );
			exit(-1);
		}
		if( fread( double_rngs, 1, sizeof( dsfmt_t ), f_input ) < sizeof( dsfmt_t ) ) {
			fprintf( stderr, "Error on reading the random real number generator state from file %s.\n", str_f_input );
			exit(-1);
		}
		fclose( f_input );
	}
	
	#pragma omp parallel for
	for( i = 0; i < nb_neurons; ++i ) {
		/* The b function */
		lambda[ i ] = 0.1;
		a[ i ] = 1.1;
		
		spiking_times[ i ] = -1.0;
	
		reset_value[ i ] = 0.0;
	
		t_last[ i ] = 0.0;
		t_last_true[ i ] = 0.0;
		y_last[ i ] = 0.0;
		
		var[ i ] = ( fabs( lambda[ i ] ) < EPSILON )? mult_const * sqrt( sigma_squared / (2 * lambda[ i ]) ): mult_const * sqrt( sigma_squared );
		threshold[ i ] = 1.0;
		max[ i ] = 0.0;
		max_cum_sum[ i ] = 0.0;
	}
fprintf( stdout, "After omp loop\n" );
fflush( stdout );
	for( i = 0; i < indices_size; ++i ) {
		for( j = 0; j < nb_neurons; ++j ) {
			indices[ i ][ j ] = j;
		}
	}
	
	for( i = 0; i < nb_neurons; ++i ) {
		/* Draw a random integer between 0 and nb_neurons included */
		nb_couplings[ i ] = 0;
		for( j = 0; j < nb_neurons; ++j ) {
			if( dsfmt_genrand_close1_open2( double_rngs ) < conn_prob ) {
				++(nb_couplings[ i ]);
			}
		}
	}
fprintf( stdout, "Before connection type part\n" );
fflush( stdout );
	switch( conn_e ) {
		case FULL:				/* All connected, including self */
		case COMPLETE:			/* All connected, except self */
		case INDEPENDENT:		/* No connections */
			break;
		case RECONSTRUCTION:	/*  */
								for( i = 0; i < threads_size; ++i ) {
									thread_indices[ i ] = i + 1;
									pthread_create( &(threads[ i ]), NULL, init_rec_fn, (void*) &(thread_indices[ i ]) );
								}
								for( i = 0; i < threads_size; ++i ) {
									pthread_join( threads[ i ], NULL );
								}
								
			break;
		case RANDOM:			/*  */
								for( i = 0; i < nb_neurons; ++i ) {
									max_ind = size;
									check_alloc( interaction_graph[ i ] = malloc( sizeof( unsigned int ) * nb_couplings[ i ] ), "interaction_graph_i" );
									for ( j = 0; j < nb_couplings[ i ]; ++j ) {
										ind_buff = get_int_using_rng( max_ind, 0 );
										interaction_graph[ i ][ j ] = indices[ 0 ][ ind_buff ];
										indices[ 0 ][ ind_buff ] = indices[ 0 ][ max_ind ];
										indices[ 0 ][ max_ind ] = interaction_graph[ i ][ j ];
										--max_ind;
									}
								}
			break;
	}	
}

void* init_rec_fn( void * ind ) {
	unsigned int m, n, max, min, *index;
fprintf( stdout, "Bonjour de %u\n", *((unsigned int*) ind) );
fflush( stdout );
	index = (unsigned int*) ind;
	
	for( m = 0; m < nb_neurons; ++m ) {
		reconstruction_graph[ *index - 1 ][ m ] = int_rngs[ *index ];
		max_ind = size;
		min		= (*index - 1) * nb_couplings[ m ] / threads_size;
		max		= *index * nb_couplings[ m ] / threads_size;
		for( n = min; n < max; ++n ) {
			get_int_using_rng( max_ind, *index );
			--max_ind;
		}
	}
	
	return NULL;
}

inline
unsigned int get_int_using_rng( const unsigned int limit, const unsigned int ind_rng ) {
	unsigned int ret_val, u, divide;

	divide = UINT_MAX / (limit + 1);	/* {0, ..., limit}, so limit+1 possibilities */
	do {
		u = sfmt_genrand_uint64( &(int_rngs[ ind_rng ]) );
		ret_val = u / divide;
	} while( ret_val > limit );
	return ret_val;
}

void save(void) {
	f_output = fopen( str_f_output, "wb" );
	f_seeds = fopen( "../Results/Seeds/seeds.bin", "wb" );
	/* Saving the characteristic parameters of this simulation */
	fwrite( &nb_neurons, 1, sizeof( unsigned int ), f_output );
	fwrite( &conn_prob, 1, sizeof( long double ), f_output );
	fwrite( &max_time, 1, sizeof( long double ), f_output );
	fwrite( &nb_itr, 1, sizeof( unsigned int ), f_output );
	fwrite( int_rngs, 1, sizeof( sfmt_t ) * int_rngs_size, f_output );
	fwrite( double_rngs, 1, sizeof( dsfmt_t ), f_output );
	
	/* Saving the state of the seed at the end of the simulation */
	fwrite( int_rngs, 1, sizeof( sfmt_t ), f_seeds );
	fwrite( double_rngs, 1, sizeof( dsfmt_t ), f_seeds );
	
	fclose( f_output );
	fclose( f_seeds );
}

void destroy(void) {
	fclose( f_results );
	
	free( int_rngs );
	free( double_rngs );
	
	free( thread_indices );

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
	
	for( i = 0; i < indices_size; ++i ) {
		free( indices[ i ] );
	}
	free( indices );
	
	free( nb_couplings );

	switch( conn_e ) {
		case FULL:				/* All connected, including self */
		case COMPLETE:			/* All connected, including self */
		case INDEPENDENT:		/* No connections */
								interaction_graph		= NULL;
								reconstruction_graph	= NULL;
			break;
		case RECONSTRUCTION:	/*  */
								for( i = 0; i < reconstruction_graph_size; ++i ) {
									free( reconstruction_graph[ i ] );
								}
								free( reconstruction_graph );
			break;
		case RANDOM:			/*  */
								for( i = 0; i < nb_neurons; ++i ) {
									free( interaction_graph[ i ] );
								}
								free( interaction_graph );
			break;
	}
}

inline
void compute_m_i_s( const time_interval_t* const time_int ) {
	/* Initializations */
	i = 0;

	/* Algorithm */
	#pragma omp parallel for
	for( i = 0; i < nb_neurons; ++i ) {
		max[ i ] = max_prob( i, time_int );
	}
	
	max_cum_sum[ 0 ] = max[ 0 ];
	for( i = 1; i < nb_neurons; ++i ) {
		max_cum_sum[ i ] = max_cum_sum[ i - 1 ] + max[ i ];
	}
}

void simulate(void) {
	/* Declarations */
	long double		delta_t,		/* The time advance since the lower bound of the time interval */
					t,				/* A buffer for a time value */
					u[2],			/* Two uniformly generated values, used for the genration of normaly distributed numbers */
					normals[2],		/* Two normaly distributed independent numbers */
					sqrt_log;		/* A buffer for an intermediate value for the computation of the normal numbers */
	time_interval_t	time_int;		/* A time interval for the computation of the Poisson process */
	unsigned int	normals_ind;	/* The index of the normal number used for the computation of the brownian motion. @see #normals */

	/* Initializations */
	u[ 0 ] = dsfmt_genrand_open_open( double_rngs ); u[ 1 ] = dsfmt_genrand_open_open( double_rngs );
	normals[ 0 ] = sqrt( -2 * log( u[ 0 ] ) ) * cos( u[ 1 ] ); normals[ 1 ] = sqrt( -2 * log( u[ 0 ] ) ) * sin( u[ 1 ] );

	time_int.lower_bound = -time_step;
	time_int.upper_bound = 0.0;
	
	normals_ind = 0;

	/* Algorithm */
NEXT_STEP:
	/* Wait for some potentials to raise high enough */
	delta_t = 0.0;
	do {
		time_int.lower_bound = time_int.upper_bound;
		time_int.upper_bound += (time_int.lower_bound + time_step > max_time)? (max_time - time_int.lower_bound): time_step;
		if( time_int.upper_bound >= max_time || nb_accepted > nb_itr ) return;
		compute_m_i_s( &time_int );
	} while( fabs( max_cum_sum[ size ] ) < EPSILON );
NEXT_SPIKE:
	delta_t += -log( dsfmt_genrand_open_close( double_rngs ) ) / max_cum_sum[ size ];	/* Generating a time step of exponential distribution E(\sum max_i) */
	if( time_int.lower_bound + delta_t > time_int.upper_bound ) goto NEXT_STEP;			/* If next potential spike is over time bound, must change of time interval */
	spiking_neuron = which_spiking();
	t = time_int.lower_bound + delta_t - t_last[ spiking_neuron ];
	t_last[ spiking_neuron ] = time_int.lower_bound + delta_t;
	if( fabs( lambda[ spiking_neuron ] ) < EPSILON ) {
		y_last[ spiking_neuron ] += sigma * normals[ normals_ind ];
	} else {
		y_last[ spiking_neuron ] = a[ spiking_neuron ]
									+ exp( -lambda[ spiking_neuron ] * t ) * (y_last[ spiking_neuron ] - a[ spiking_neuron ])
									+ sqrt( sigma_squared * (1 - exp( -2 * lambda[ spiking_neuron ] * t )) / (2 * lambda[ spiking_neuron ]) ) * normals[ normals_ind ];
	}
fprintf( stdout, "Before drawing normal numbers\n" );
fflush( stdout );
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
fprintf( stdout, "Before rejecting spike\n" );
fflush( stdout );
	if( !(probability( spiking_neuron ) > dsfmt_genrand_close_open( double_rngs ) * max[ spiking_neuron ]) ) {
		/* The spike is rejected */
		++nb_rejected;
		goto NEXT_SPIKE;
	}
	/* From now on the jump has been accepted,  */
	++ nb_accepted;
	y_last[ spiking_neuron ] = reset_value[ spiking_neuron ];

	spiking_times[ spiking_times_array_index ] = time_int.lower_bound + delta_t;
	++spiking_times_array_index;
	if( spiking_times_array_index > size ) {
		fwrite( spiking_times, sizeof( long double ), nb_neurons, f_results );
		spiking_times_array_index = 0;
	}
fprintf( stdout, "Before reconstruction\n" );
fflush( stdout );
	switch( conn_e ) {
		case RANDOM:			/* The indices of the postsynaptic neurons are stored in the interaction_graph array */
								#pragma omp parallel for
								for( j = 0; j < nb_couplings[ spiking_neuron ]; ++j ) {
									y_last[ interaction_graph[ spiking_neuron ][ j ] ] += interaction( spiking_neuron, j );
								}
			break;
		case FULL:				/* All neurons are postsynaptic neurons */
								#pragma omp parallel for
								for( j = 0; j < nb_neurons; ++j ) {
										y_last[ j ] += interaction( spiking_neuron, j );
								}
			break;
		case COMPLETE:			/* All neurons are postsynaptic neurons, except spiking neuron */
								#pragma omp parallel for
								for( j = 0; j < nb_neurons; ++j ) {
									if( j != spiking_neuron ) {
										y_last[ j ] += interaction( spiking_neuron, j );
									}
								}
			break;
		case INDEPENDENT:		/* There are no connections between neurons, nothing to be done in this case. */
			break;
		case RECONSTRUCTION:	/* Interaction graph was too big, must be reconstructed on the fly */
								for( i = 0; i < threads_size; ++i ) {
									thread_indices[ i ] = i + 1;
									pthread_create( &(threads[ i ]), NULL, reconstruction_fn, (void*) &(thread_indices[ i ]) );
								}
								for( i = 0; i < threads_size; ++i ) {
									pthread_join( threads[ i ], NULL );
								}
								
								/* Sorting the array of indices after use */
								#pragma omp parallel for
								for( i = 0; i < NB_THREADS; ++i ) {
									for( j = 0; j < nb_neurons; ++j ) {
										indices[ i ][ j ] = j;
									}
								}
			break;
	}

	time_int.upper_bound = time_int.lower_bound + delta_t;
	goto NEXT_STEP;
}

void* reconstruction_fn( void * ind ) {
	unsigned int m, n, max, min, *index;
	
	index = (unsigned int*) ind;
	
	int_rngs[ *index ] = reconstruction_graph[ *index - 1 ][ spiking_neuron ];
	/* At this point it is assumed the array of indices 'indices' is sorted */
	max_ind	= size;
	min		= (*index - 1) * nb_couplings[ spiking_neuron ] / threads_size;
	max		= *index * nb_couplings[ spiking_neuron ] / threads_size;
	for( m = min; m < max; ++m ) {
		n = get_int_using_rng( max_ind, *index ); /* A random non already chosen index in the array of indices (meta) */
		y_last[ indices[ *index - 1 ][ n ] ] += interaction( spiking_neuron, n );

		/* The chosen index is put at the end of the array so that it will not be chosen again */
		ind_buff = indices[ *index - 1 ][ n ];
		indices[ *index - 1 ][ n ] = indices[ *index ][ max_ind ];
		indices[ *index - 1 ][ max_ind ] = ind_buff;
		/* The maximum index of the never-chosen indices is decreased as there is a new chosen index */
		--max_ind;
	}
	
	return NULL;
}

inline
long double interaction( const unsigned int presynaptic, const unsigned int postsynaptic ) {
	(void)presynaptic; (void)postsynaptic; /* The parameters are not used currently */
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
	if( fabs( lambda[ neuron ] ) < EPSILON ) {
		max_prob = y_last[ neuron ] + var[ neuron ]
				 - threshold[ neuron ];
	} else {
		max_prob = a[ neuron ]
					+ exp( -lambda[ neuron ] * (time_int->upper_bound - t_last[ neuron ]) ) * (y_last[ neuron ] - a[ neuron ])
					+ var[ neuron ] * sqrt(1 - exp( -2 * lambda[neuron] * (time_int->upper_bound - t_last[neuron]) ))
				 - threshold[ neuron ];
	}
	return (max_prob > 0)? slope * pow( max_prob, exponent ): 0.0;
}

