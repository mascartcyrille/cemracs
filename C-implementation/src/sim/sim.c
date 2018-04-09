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
#include <tgmath.h>
#include <limits.h>
#include <time.h>
#include <sys/stat.h>
#include <omp.h>
#include <string.h>

#include "sim.h"
#include "../lib/SFMT-src-1.5.1/SFMT.h"		/* Simple precision rng */
#include "../lib/dSFMT-src-2.2.3/dSFMT.h"	/* Double precision rng */

/* DEFINES */
#define 	MAX_MEMORY	1e0					/* The maximum memory allowed for the storing of the interaction graph. */
#define 	BOOL		unsigned char		/* Boolean type, "equivalent" to the one defined in ++c */
#define 	EPSILON		1e-15				/* The minimum comparison value for floating point numbers.
										   		Any real number of absolute value below EPSILON is considered equal to zero.
											 */
#define		MAXSIZE		256					/* The maximum number of characters */

/* GLOBAL VARIABLES */
time_t			rawtime, init_time;				/* Raw time info (in seconds, since 01/01/1970) sample at the beginning of the simulation */
struct tm	*	timeinfo;					/* The time info in human redable form */

unsigned int	i, j,						/* Two generic buffer indices */
				ind_buff, max_ind,			/* Indices for the array of indices (meta) */
				spiking_times_array_index,	/* The current memory box available for storing a time value */
				nb_neurons,					/* Denoted N hereafter */
				size,						/* = nb_neurons - 1 */
				nb_itr,						/* Maximum number of jumps (sampling acceptations) allowed for the simulation */
				nb_accepted,				/* Number of accepted jumps */
				nb_rejected;				/* Number of rejected jumps */

unsigned char	is_eof,						/* Is it the end of the line (text f_params) */
				is_eol;						/* Is it the end of the file (text f_params) */

long double	conn_prob,	/* Probability of a connection between two neurons, possibly self */
			max_time,	/* The maximum allowed time. */
/* Probability of spike function: slope * (Pos(potential-threshold))^exponent */
			slope,		/* The slope as defined above */
			exponent,	/* The exponent as defined above */
			time_step,	/* For the thinning method, the size of the time intervals in which spikes are looked for */
/* Max variance for Brownian motion */
			sigma,			/* A control for the importance of the Brownian motion */
			confidence,		/* A constant, caculated from, but not only, sigma */
			inter,			/* The weight of the interaction */
			beta;			/*  */

conn_type	conn_e;		/* What is the type of connection of this simulation */
param_file	param_e;	/* The encoding for the parameter file */

/* The rngs and associated threads */
sfmt_t		*	int_rngs,			/* Array [1..NB_THREADS]: one per thread plus one for the main thread */
			*	init_int_rngs;		/* The init value of int_rngs */
dsfmt_t		*	double_rngs,		/* Array [1..NB_THREADS]: one per thread plus one for the main thread */
			*	init_double_rngs;	/* The init value of double_rngs */
uint32_t	*	int_seed,			/* A seed buffer for the integer rng initialisation */
			*	double_seed;		/* A seed buffer for the floating point number rng initialisation */

/* Seeds, parameters and results files */
FILE		*	f_params,	/*  */
			*	f_results,	/*  */
			*	f_stats,	/*  */
			*	f_seeds,	/*  */
			*	f_reproc;	/*  */
fpos_t			eol;		/* The position of the end of the line (for use while reading f_params, both in txt and bin file) */

/* The b function */
long double 	*	lambda,	/* Array [1..N]: lambda constant values */
				*	a,		/* Array [1..N]: constant value at origin for the b function (see header) */
/* The spiking times */
				*	spiking_times,	/* Array [1..N]: the times of true spikes. The array is written down to the result file when full */
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
char		c,						/* A buffer for storing character values while reading a text file */
		*	str_f_params,			/* The name of the file of parameters */
			str_f_results[MAXSIZE],	/* The name of the file of  */
			str_f_stats[MAXSIZE],	/* The name of the file for statistics */
			str_f_seeds[MAXSIZE],	/*  */
			str_f_reproc[MAXSIZE],	/* The name of the file of  */
			str_folder[MAXSIZE],	/* The name of the folder containing all these files */
			buff[MAXSIZE];			/* A buffer of words, used during the parsing of the parameters */

/* FUNCTIONS IMPLEMENTATION */
int main( int argc, char** const argv ) {
	time( &rawtime );
	timeinfo = localtime( &rawtime );

	/* Sets if binary or text parameter file */
	if( argv[ 1 ][ 0 ] == '-' && argv[ 1 ][ 1 ] == 'b' ) {
		param_e = BIN;
	} else if( argv[ 1 ][ 0 ] == '-' && argv[ 1 ][ 1 ] == 't' ) {
		param_e = TXT;
	} else {
		param_e = UNKNOWN;
		fprintf( stderr, "\nError while trying to detect if the parameter file is of type text of binary. Aborting...\n\n" );
		exit( -1 );
	}
	/* Checks the number of arguments */
	switch( argc ) { 
		case 3:		/* A file that contains all parameters, text encoded if it is new or binary encoded if it is for reproduction */
					str_f_params = argv[ 2 ];
					parse_base();					/* Gets the number of neurons and probability of connection for this simulation. */
					create();						/* Dynamically allocates memory. Should be called only once, except if the number of neurons changes. */
					init();							/* Initializes the variables values. */
					time( &init_time );
					fprintf( stdout, "Initialization time: %u\n", (unsigned int) ( init_time - rawtime ) );
					simulate();						/* Makes the simulation. */
					if( param_e == TXT ) save();	/* Saves the parameters in files, so that the simulation can be redone. Not called if this is actually a redoing simulation. */
					destroy();						/* Frees the memory dynamically allocated. */
			break;
		default:	/* Either too few or too many arguments. In either cases, send an error message and ends the simulation. */
					fprintf( stderr, "Bad arguments.\n\tsim -[t|b] param_file\n" );
					exit( -1 );
			break;
	}

	return 0;
}

inline
void parse_base() {
	unsigned char found;
	switch( param_e ) {
		case TXT:	/* If the parameters are in a text file */
					f_params = fopen( str_f_params, "r" );
					if( f_params != NULL ) {
						found = 0;
						do {
							/* Get key */
							get_value( '=' );
							/* Parse */
							if( strcmp( buff, "nb_neurons" ) == 0 ) {
								get_value( ',' );
								nb_neurons = strtoumax( buff, NULL, 10 );
								++found;
							} else if( strcmp( buff, "conn_prob" ) == 0 ) {
								get_value( ',' );
								conn_prob = strtold( buff, NULL );
								++found;
							}
							if( errno != 0 ) {
								fprintf( stderr, "%s\n", strerror( errno ) );
								exit( -1 );
							}
						} while( c != EOF && found != 2 );
						fclose( f_params );
					} else {
						fprintf( stderr, "\nCannot open the text file %s of parameters.\n\n", str_f_params );
						exit( -1 );
					}
			break;
		case BIN:	/* The parameter file is binary */
					f_params = fopen( str_f_params, "rb" );
					if( f_params != NULL ) {
						if( fread( &nb_neurons, sizeof( unsigned int ), 1, f_params ) < 1 ) {
							fprintf( stderr, "\nCannot read the number of neurons in binary parameter file %s. Aborting...\n\n", str_f_params );
							exit( -1 );
						}
						if( fread( &conn_prob, sizeof( long double ), 1, f_params ) < 1 ) {
							fprintf( stderr, "\nCannot read the probability of connection in binary parameter file %s. Aborting...\n\n", str_f_params );
							exit( -1 );
						}
						fgetpos( f_params, &eol );
						fclose( f_params );
					} else {
						fprintf( stderr, "\nCannot open the text file %s of parameters.\n\n", str_f_params );
						exit( -1 );
					}
			break;
		default:	/* Shall never enter here */
					fprintf( stderr, "\nUnknown error: entered in the default case in parse_base function. Aborting...\n\n" );
					exit( -1 );
			break;
	}
}

inline
void get_value( const char sep ) {
	i = 0; c = fgetc( f_params );
	while( i < (MAXSIZE - 1) && c != EOF && c != '\n' && c != sep ) {
		if( c != ' ' && c != '\t' ) {
			buff[i] = c;
			++i;
		}
		c = fgetc( f_params );
	}
	is_eol = ( c == '\n' );
	is_eof = ( c == EOF );
	buff[i] = '\0';
}

inline
void set_value_ld( long double * value ) {
	fpos_t			bol;
	unsigned int	line_size,
					nb_words;

	line_size = 1;
	nb_words = 0;

	fgetpos( f_params, &bol );	/* Get the beginning of line position */
	/* Read until the end of the line */
	for( c = fgetc( f_params ); c != EOF && c != '\n'; c = fgetc( f_params ) ) {
		line_size += ( c == ',' )? 1: 0;
	}
	fgetpos( f_params, &eol );	/* Get the end of line position */
	fsetpos( f_params, &bol );	/* Set the position to the beginning of the line */
	do {
		get_value( ',' );	/* Get the value */
		if( nb_words % line_size == 0 ) {	/* Reset cursor to the beginning of the line if all values have been used */
			fsetpos( f_params, &bol );
		}
		value[ nb_words++ ] = strtold( buff, NULL );	/* Sets the value for a */
	} while( nb_words < nb_neurons && !is_eof );
	fsetpos( f_params, &eol );	/* Sets the position in the file at the end of the line, so that the parsing can continue */
}

inline
void parse( char * const str_file ) {
	size_t fread_size = 0;

	str_f_params = str_file;
	switch( param_e ) {
		case TXT:	/* Text parameter file */
					f_params = fopen( str_f_params, "r" );
					if( f_params != NULL ) {
						do {
							/* Get key */
							get_value( '=' );
							/* Parse */
							if( strcmp( buff, "nb_itr" ) == 0 ) {
								get_value( ',' );
								nb_itr = strtoumax( buff, NULL, 10 );
							} else if( strcmp( buff, "max_time" ) == 0 ) {
								get_value( ',' );
								max_time = strtold( buff, NULL );
							} else if( strcmp( buff, "slope" ) == 0 ) {
								get_value( ',' );
								slope = strtold( buff, NULL );
							} else if( strcmp( buff, "exponent" ) == 0 ) {
								get_value( ',' );
								exponent = strtold( buff, NULL );
							} else if( strcmp( buff, "time_step" ) == 0 ) {
								get_value( ',' );
								time_step = strtold( buff, NULL );
							} else if( strcmp( buff, "sigma" ) == 0 ) {
								get_value( ',' );
								sigma = strtold( buff, NULL );
							} else if( strcmp( buff, "confidence" ) == 0 ) {
								get_value( ',' );
								confidence = strtold( buff, NULL );
							} else
								if( strcmp( buff, "a" ) == 0 ) {
								set_value_ld( a );
							} else if( strcmp( buff, "lambda" ) == 0 ) {
								set_value_ld( lambda );
							} else if( strcmp( buff, "reset_value" ) == 0 ) {
								set_value_ld( reset_value );
							} else if( strcmp( buff, "threshold" ) == 0 ) {
								set_value_ld( threshold );
							} else
								if( strcmp( buff, "f_reproc" ) == 0 ) {
								get_value( '\0' );
								strcpy( str_f_reproc, buff );
							} else if( strcmp( buff, "f_results" ) == 0 ) {
								get_value( '\0' );
								strcpy( str_f_results, buff );
							} else if( strcmp( buff, "f_stats" ) == 0 ) {
								get_value( '\0' );
								strcpy( str_f_stats, buff );
							} else if( strcmp( buff, "f_seeds" ) == 0 ) {
								get_value( '\0' );
								strcpy( str_f_seeds, buff );
							} else if( strcmp( buff, "folder" ) == 0 ) {
								get_value( '\0' );
								strcpy( str_folder, buff );
							}
						} while( !is_eof );
						fclose( f_params );
					} else {
						fprintf( stderr, "\nError while opening the parameter file %s.\n\n", str_f_params );
						exit( -1 );
					}
			break;
		case BIN:	/* The parameters are in a binary file */
					f_params = fopen( str_f_params, "rb" );
					if( f_params != NULL ) {
						/* Series of fread, must be in the same order than the fwrite of the save() function */
						fsetpos( f_params, &eol );
						fread_size = fread( &nb_itr, sizeof( unsigned int ), 1, f_params );
						if( fread_size < 1 ) {
							fprintf( stderr, "\nError while trying to get parameter nb_itr from file of parameters %s. Aborting...\n\n", str_f_params );
							exit( -1 );
						}

						fread_size = fread( &max_time, sizeof( long double ), 1, f_params );
						if( fread_size < 1 ) {
							fprintf( stderr, "\nError while trying to get parameter max_time from file of parameters %s. Aborting...\n\n", str_f_params );
							exit( -1 );
						}
						fread_size = fread( &slope, sizeof( long double ), 1, f_params );
						if( fread_size < 1 ) {
							fprintf( stderr, "\nError while trying to get parameter slope from file of parameters %s. Aborting...\n\n", str_f_params );
							exit( -1 );
						}
						fread_size = fread( &exponent, sizeof( long double ), 1, f_params );
						if( fread_size < 1 ) {
							fprintf( stderr, "\nError while trying to get parameter exponent from file of parameters %s. Aborting...\n\n", str_f_params );
							exit( -1 );
						}
						fread_size = fread( &time_step, sizeof( long double ), 1, f_params );
						if( fread_size < 1 ) {
							fprintf( stderr, "\nError while trying to get parameter time_step from file of parameters %s. Aborting...\n\n", str_f_params );
							exit( -1 );
						}
						fread_size = fread( &sigma, sizeof( long double ), 1, f_params );
						if( fread_size < 1 ) {
							fprintf( stderr, "\nError while trying to get parameter sigma from file of parameters %s. Aborting...\n\n", str_f_params );
							exit( -1 );
						}
						fread_size = fread( &confidence, sizeof( long double ), 1, f_params );
						if( fread_size < 1 ) {
							fprintf( stderr, "\nError while trying to get parameter confidence from file of parameters %s. Aborting...\n\n", str_f_params );
							exit( -1 );
						}
						
						fread_size = fread( a, sizeof( long double ), nb_neurons, f_params );
						if( fread_size < nb_neurons ) {
							fprintf( stderr, "\nError while trying to get parameter \"a\" from file of parameters %s. Aborting...\n\n", str_f_params );
							exit( -1 );
						}
						fread_size = fread( lambda, sizeof( long double ), nb_neurons, f_params );
						if( fread_size < nb_neurons ) {
							fprintf( stderr, "\nError while trying to get parameter lambda from file of parameters %s. Aborting...\n\n", str_f_params );
							exit( -1 );
						}
						fread_size = fread( reset_value, sizeof( long double ), nb_neurons, f_params );
						if( fread_size < nb_neurons ) {
							fprintf( stderr, "\nError while trying to get parameter reset_value from file of parameters %s. Aborting...\n\n", str_f_params );
							exit( -1 );
						}
						fread_size = fread( threshold, sizeof( long double ), nb_neurons, f_params );
						if( fread_size < nb_neurons ) {
							fprintf( stderr, "\nError while trying to get parameter threshold from file of parameters %s. Aborting...\n\n", str_f_params );
							exit( -1 );
						}

						fread_size = fread( str_f_reproc, sizeof( char ), MAXSIZE, f_params );
						if( fread_size < MAXSIZE ) {
							fprintf( stderr, "\nError while trying to get parameter str_f_reproc from file of parameters %s. Aborting...\n\n", str_f_params );
							exit( -1 );
						}
						fread_size = fread( str_f_results, sizeof( char ), MAXSIZE, f_params );
						if( fread_size < MAXSIZE ) {
							fprintf( stderr, "\nError while trying to get parameter str_f_results from file of parameters %s. Aborting...\n\n", str_f_params );
							exit( -1 );
						}
						fread_size = fread( str_f_stats, sizeof( char ), MAXSIZE, f_params );
						if( fread_size < MAXSIZE ) {
							fprintf( stderr, "\nError while trying to get parameter str_f_stats from file of parameters %s. Aborting...\n\n", str_f_params );
							exit( -1 );
						}
						fread_size = fread( str_f_seeds, sizeof( char ), MAXSIZE, f_params );
						if( fread_size < MAXSIZE ) {
							fprintf( stderr, "\nError while trying to get parameter str_f_seeds from file of parameters %s. Aborting...\n\n", str_f_params );
							exit( -1 );
						}
						fread_size = fread( str_folder, sizeof( char ), MAXSIZE, f_params );
						if( fread_size < MAXSIZE ) {
							fprintf( stderr, "\nError while trying to get parameter str_folder from file of parameters %s. Aborting...\n\n", str_f_params );
							exit( -1 );
						}

						fgetpos( f_params, &eol );
						fclose( f_params );
					} else {
						fprintf( stderr, "\nError while opening the parameter file %s.\n\n", str_f_params );
						exit( -1 );
					}
			break;
		default:	/* Should never enter here */
					fprintf( stderr, "\nError while .\n\n" );
					exit( -1 );
			break;
	}
}

inline
void check_alloc( const void* const var, const char * const var_name ) {
	if( var == NULL ) {
		fprintf( stderr, "Cannot allocate variable %s\n", var_name );
		abort();
	}
}

void create(void) {
	conn_e =	( conn_prob >= (long double) 1.0 )?									FULL:
				( conn_prob <= (long double) 0.0 )?									INDEPENDENT:
				( conn_prob * nb_neurons * nb_neurons > (long double) MAX_MEMORY )?	RECONSTRUCTION:
																					RANDOM;

	check_alloc( int_rngs			= malloc( sizeof( sfmt_t ) ), "integer rngs" );
	check_alloc( init_int_rngs		= malloc( sizeof( sfmt_t ) ), "initial integer rngs" );
	check_alloc( double_rngs		= malloc( sizeof( dsfmt_t ) ), "double rngs" );
	check_alloc( init_double_rngs	= malloc( sizeof( dsfmt_t ) ), "initial double rngs" );
	
	check_alloc( lambda				= malloc( sizeof( long double ) * nb_neurons ), "lambda" );
	check_alloc( a					= malloc( sizeof( long double ) * nb_neurons ), "a" );
	
	check_alloc( spiking_times		= malloc( sizeof( long double ) * nb_neurons ), "spiking_times" );
	check_alloc( reset_value		= malloc( sizeof( long double ) * nb_neurons ), "reset_value" );
	check_alloc( t_last				= malloc( sizeof( long double ) * nb_neurons ), "t_last" );
	check_alloc( t_last_true		= malloc( sizeof( long double ) * nb_neurons ), "t_last_true" );
	check_alloc( y_last				= malloc( sizeof( long double ) * nb_neurons ), "y_last" );
	check_alloc( var				= malloc( sizeof( long double ) * nb_neurons ), "var" );
	check_alloc( threshold			= malloc( sizeof( long double ) * nb_neurons ), "threshold" );
	check_alloc( max				= malloc( sizeof( long double ) * nb_neurons ), "max" );
	check_alloc( max_cum_sum		= malloc( sizeof( long double ) * nb_neurons ), "max_cum_sum" );
	check_alloc( indices			= malloc( sizeof( long double ) * nb_neurons ), "indices" );
	
	check_alloc( nb_couplings = malloc( sizeof( unsigned int ) * nb_neurons ), "nb_couplings" );
	
	switch( conn_e ) {
		case FULL:				/* Graph fully connected, including self */
		case COMPLETE:			/* Complete graph (all connected to all, except self) */
		case INDEPENDENT:		/* No connections */
								interaction_graph		= NULL;
								reconstruction_graph	= NULL;
			break;
		case RECONSTRUCTION:	/* The graph will be reconstructed on-the-fly during simulation */
								interaction_graph		= NULL;
								check_alloc( reconstruction_graph = malloc( sizeof( sfmt_t ) * nb_neurons ), "reconstruction_graph" );
			break;
		case RANDOM:			/* The graph is "random" and small enough that it can be fully stored in the RAM */
								reconstruction_graph	= NULL;
								check_alloc( interaction_graph = malloc( sizeof( unsigned int* ) * nb_neurons ), "interaction_graph" );
			break;
	}
}

void files_and_folders( void ) {
	sprintf( buff, "./results/params/%s", str_f_seeds );	sprintf( str_f_seeds, "%s", buff );
	char * res_folder = "./results";

	strcpy( buff, str_folder );
	sprintf( str_folder, "%s/%s-%04d-%02d-%02d-%02d-%02d-%02d", res_folder, buff, timeinfo->tm_year + 1900, timeinfo->tm_mon + 1, timeinfo->tm_mday, timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec );
	i = 0;
	do {
		do {
			buff[ i ] = str_folder[ i ];
			++i;
		} while( str_folder[ i ] != '/' && str_folder[ i ] != '\0' );
		buff[ i ] = '\0';
		mkdir( buff, 0700 );
		buff[ i ] = str_folder[ i ];
		++i;
	} while( str_folder[ i ] != '\0' );

	sprintf( buff, "%s/%s", str_folder, str_f_reproc );		sprintf( str_f_reproc, "%s", buff );
	sprintf( buff, "%s/%s", str_folder, str_f_stats );		sprintf( str_f_stats, "%s", buff );
	sprintf( buff, "%s/%s", str_folder, str_f_results );	sprintf( str_f_results, "%s", buff );
}

void init(void) {
	parse( str_f_params );	/* Initializes the variables of arbitraty values. */
	
	size			= nb_neurons - 1;	/* For manipulating the arrays of size nb_neurons. */
	nb_accepted		= 0;				/* Number of accepted points during the thinning method. */
	nb_rejected		= 0;				/* Number of rejected points during the thinning method. */
	spiking_times_array_index	= 0;	/* Pointer to the last used memory box of the spiking_times array */

	/* Seeds the RNGs */
	switch( param_e ) {
		case TXT:	/* Parameters from a text file */
					files_and_folders();	/* Initializes the file names. */
					f_seeds = fopen( str_f_seeds, "rb" );
					if( f_seeds != NULL ) {
						if( fread( init_int_rngs, 1, sizeof( sfmt_t ), f_seeds ) < sizeof( sfmt_t ) ) {
							fprintf( stderr, "Error on reading the random integer generator state from file %s.\n", str_f_seeds );
							exit(-1);
						} else {
							*int_rngs = *init_int_rngs;
						}
						if( fread( init_double_rngs, 1, sizeof( dsfmt_t ), f_seeds ) < sizeof( dsfmt_t ) ) {
							fprintf( stderr, "Error on reading the random real number generator state from file %s.\n", str_f_seeds );
							exit(-1);
						} else {
							*double_rngs = *init_double_rngs;
						}
						fclose( f_seeds );
					} else {
						fprintf( stderr, "Error while opening seed file. Using current time value for initialisation.\n" );
						sfmt_init_gen_rand( init_int_rngs, rawtime );
						*int_rngs = *init_int_rngs;
						dsfmt_init_gen_rand( init_double_rngs, rawtime );
						*double_rngs = *init_double_rngs;
					}
			break;
		case BIN:	/* Parameters from a binary file */
					f_params = fopen( str_f_params, "rb" );
					if( f_params != NULL ) {
						fsetpos( f_params, &eol );
						if( fread( init_int_rngs, 1, sizeof( sfmt_t ), f_params ) < sizeof( sfmt_t ) ) {
							fprintf( stderr, "Error on reading the random integer generator state from file %s.\n", str_f_params );
							exit(-1);
						}
						if( fread( init_double_rngs, 1, sizeof( dsfmt_t ), f_params ) < sizeof( dsfmt_t ) ) {
							fprintf( stderr, "Error on reading the random real number generator state from file %s.\n", str_f_params );
							exit(-1);
						}
						*int_rngs = *init_int_rngs;
						*double_rngs = *init_double_rngs;
						fclose( f_params );
					} else {
						fprintf( stderr, "\nError while reading the seed in reproducibility file %s. Aborting...\n\n", str_f_params );
						exit( -1 );
					}
			break;
		default:	/* Should never enter here */
					fprintf( stderr, "\nDefault case of init switch activated. Aborting...\n\n" );
					exit( -1 );
			break;
	}
	
	/* The b function */
	for( i = 0; i < nb_neurons; ++i ) {
		spiking_times[ i ]	=	-1.0;

		t_last[ i ]			=	0.0;
		t_last_true[ i ]	=	0.0;
		y_last[ i ]			=	0.0;

		var[ i ]			=	( fabsl( lambda[ i ] ) < EPSILON )?	confidence * sigma * 1 / sqrtl( 2 * lambda[ i ] ):
																	confidence * sigma;
		max[ i ]			=	0.0;
		max_cum_sum[ i ]	=	0.0;

		indices[ i ]		=	i;
	}
	
	if( conn_e == RECONSTRUCTION || conn_e == RANDOM ) {	/* Generates the couplings only of the graph will be used */
		for( i = 0; i < nb_neurons; ++i ) {
			/* Draw a random integer between 0 and nb_neurons included */
			nb_couplings[ i ] = 0;
			for( j = 0; j < nb_neurons; ++j ) {
				if( dsfmt_genrand_close1_open2( double_rngs ) < conn_prob ) {
					++(nb_couplings[ i ]);
				}
			}
		}
	}
	
	switch( conn_e ) {
		case FULL:				/* All connected, including self */
		case COMPLETE:			/* All connected, except self */
		case INDEPENDENT:		/* No connections */
								fprintf( stdout, "Chosen the FULL, COMPLETE or INDEPENDENT connection type.\n" );
			break;
		case RECONSTRUCTION:	/* The graph is too big and must be reconstructed on-the-fly */
								fprintf( stdout, "Chosen the RECONSTRUCTION connection type.\n" );
								for( i = 0; i < nb_neurons; ++i ) {
									reconstruction_graph[ i ] = *int_rngs;
									max_ind = size;
									for ( j = 0; j < nb_couplings[ i ]; ++j ) {
										get_int( max_ind );
										--max_ind;
									}
								}
			break;
		case RANDOM:			/* The graph is small enough to be stored in memory */
								fprintf( stdout, "Chosen the RANDOM connection type.\n" );
								for( i = 0; i < nb_neurons; ++i ) {
									max_ind = size;
									if( interaction_graph[ i ] != NULL ) {
										free( interaction_graph[ i ] );
									}
									check_alloc( interaction_graph[ i ] = malloc( sizeof( unsigned int ) * nb_couplings[ i ] ), "interaction_graph_i" );
									for ( j = 0; j < nb_couplings[ i ]; ++j ) {
										ind_buff = get_int( max_ind );
										interaction_graph[ i ][ j ] = indices[ ind_buff ];
										indices[ ind_buff ] = indices[ max_ind ];
										indices[ max_ind ] = interaction_graph[ i ][ j ];
										--max_ind;
									}
								}
			break;
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
	/* Saving some statistics about the simulation (accepted/rejected number of jumps, etc.) */
	f_stats = fopen( str_f_stats, "w" );
	if( f_stats != NULL ) {
		fprintf( f_stats, "Accepted spikes: %d\n", nb_accepted );
		fprintf( f_stats, "Rejected spikes: %d\n", nb_rejected );
		fprintf( f_stats, "Number of spikes: %d\n", nb_itr );
		fclose( f_stats );
	} else {
		fprintf( stderr, "\nError while writing down the statistics of the simulation in the statistic file %s.\n\n", str_f_stats );
	}

	/* Saving the characteristic parameters of this simulation */
	f_reproc = fopen( str_f_reproc, "wb" );
	if( f_reproc != NULL ) {
		fwrite( &nb_neurons, sizeof( unsigned int ), 1, f_reproc );
		fwrite( &conn_prob, sizeof( long double ), 1, f_reproc );

		fwrite( &nb_itr, sizeof( unsigned int ), 1, f_reproc );
		fwrite( &max_time, sizeof( long double ), 1, f_reproc );
		fwrite( &slope, sizeof( long double ), 1, f_reproc );
		fwrite( &exponent, sizeof( long double ), 1, f_reproc );
		fwrite( &time_step, sizeof( long double ), 1, f_reproc );
		fwrite( &sigma, sizeof( long double ), 1, f_reproc );
		fwrite( &confidence, sizeof( long double ), 1, f_reproc );

		fwrite( a, sizeof( long double ), nb_neurons, f_reproc );
		fwrite( lambda, sizeof( long double ), nb_neurons, f_reproc );
		fwrite( reset_value, sizeof( long double ), nb_neurons, f_reproc );
		fwrite( threshold, sizeof( long double ), nb_neurons, f_reproc );

		fwrite( str_f_reproc, sizeof( char ), MAXSIZE, f_reproc );
		fwrite( str_f_results, sizeof( char ), MAXSIZE, f_reproc );
		fwrite( str_f_stats, sizeof( char ), MAXSIZE, f_reproc );
		fwrite( str_f_seeds, sizeof( char ), MAXSIZE, f_reproc );
		fwrite( str_folder, sizeof( char ), MAXSIZE, f_reproc );

		fwrite( init_int_rngs, sizeof( sfmt_t ), 1, f_reproc );
		fwrite( init_double_rngs, sizeof( dsfmt_t ), 1, f_reproc );
		fclose( f_reproc );
	} else {
		fprintf( stderr, "\nError while writing down the parameters of the simulation in the reproducibility file %s.\n\n", str_f_reproc );
	}

	/* Saving the state of the seed at the end of the simulation */
	mkdir( "./results/params/", 0700 );
	f_seeds = fopen( "./results/params/seeds.bin", "wb" );
	if( f_seeds != NULL ) {
		fwrite( int_rngs, 1, sizeof( sfmt_t ), f_seeds );
		fwrite( double_rngs, 1, sizeof( dsfmt_t ), f_seeds );
		fclose( f_seeds );
	} else {
		fprintf( stderr, "\nError while writing down the seeds values of the rngs in the seed file %s.\n\n", str_f_seeds );
	}
}

void destroy(void) {
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

	switch( conn_e ) {
		case FULL:				/* Fully connected graph, including self */
		case COMPLETE:			/* Complete graph (no self-connection) */
		case INDEPENDENT:		/* No connection at all */
								interaction_graph		= NULL;
								reconstruction_graph	= NULL;
			break;
		case RECONSTRUCTION:	/* Graph too big, reconstructed on-the-fly */
								free(reconstruction_graph);
			break;
		case RANDOM:			/* In memory "random" graph */
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
	sfmt_t			save_rng;		/* An RNG state, used before the reconstruction of the interaction graph */
	long double		delta_t,		/* The time advance since the lower bound of the time interval */
					tmp,			/*  */
					t,				/* A buffer for a time value */
					u[2],			/* Two uniformly generated values, used for the genration of normaly distributed numbers */
					normals[2],		/* Two normaly distributed independent numbers */
					sqrt_log;		/* A buffer for an intermediate value for the computation of the normal numbers */
	time_interval_t	time_int;		/* A time interval for the computation of the Poisson process */
	unsigned int	spiking_neuron,	/* The index of the current spiking neuron */
					normals_ind,	/* The index of the normal number used for the computation of the brownian motion. @see #normals */
					nb_itr_per_cent;

	/* Initializations */
	f_results		= fopen( str_f_results, "wb" );
	u[ 0 ]			= dsfmt_genrand_open_open( double_rngs ); u[ 1 ] = dsfmt_genrand_open_open( double_rngs );
	normals[ 0 ]	= sqrtl( -2 * logl( u[ 0 ] ) ) * cos( u[ 1 ] ); normals[ 1 ] = sqrtl( -2 * logl( u[ 0 ] ) ) * sin( u[ 1 ] );

	time_int.lower_bound = -time_step;
	time_int.upper_bound = 0.0;
	
	normals_ind		= 0;
	nb_itr_per_cent	= nb_itr / 100;
	/* Algorithm */
NEXT_STEP:
	/* Wait for some potentials to raise high enough */
	delta_t = 0.0;
	do {
		time_int.lower_bound = time_int.upper_bound;
		time_int.upper_bound += (time_int.lower_bound + time_step > max_time)? (max_time - time_int.lower_bound): time_step;
		if( time_int.upper_bound >= max_time || nb_accepted >= nb_itr ) goto END;
		compute_m_i_s( &time_int );
	} while( fabsl( max_cum_sum[ size ] ) < EPSILON );
NEXT_SPIKE:
	tmp = -logl( dsfmt_genrand_open_close( double_rngs ) ) / max_cum_sum[ size ];	/* Generating a time step of exponential distribution E(\sum max_i) */
	if( tmp <= (long double) 0.0 ) {
		fprintf( stderr, "Time advance of zero!\n" );
	}
	delta_t += tmp;
	if( time_int.lower_bound + delta_t > time_int.upper_bound ) goto NEXT_STEP;		/* If next potential spike is over time bound, must change of time interval */
	spiking_neuron = which_spiking();
	t = time_int.lower_bound + delta_t - t_last[ spiking_neuron ];
	t_last[ spiking_neuron ]	= time_int.lower_bound + delta_t;
	if( fabsl( lambda[ spiking_neuron ] ) < EPSILON ) {
		y_last[ spiking_neuron ]   += sigma * normals[ normals_ind ];
	} else {
		y_last[ spiking_neuron ]	= a[ spiking_neuron ]
									+ expl( -lambda[ spiking_neuron ] * t ) * (y_last[ spiking_neuron ] - a[ spiking_neuron ])
									+ sigma * sqrtl( (1 - expl( -2 * lambda[ spiking_neuron ] * t )) / (2 * lambda[ spiking_neuron ]) ) * normals[ normals_ind ];
	}
	if( normals_ind ) {
		/* Marsaglia's method for generating normally distributed independent random numbers */
		do {
			u[ 0 ]		= 2.0 * dsfmt_genrand_open_open( double_rngs ) - 1.0;	/* cos( U([0,1]) ) ~ U([-1,1]) */
			u[ 1 ]		= 2.0 * dsfmt_genrand_open_open( double_rngs ) - 1.0;	/* sin( U([0,1]) ) ~ U([-1,1]) */
			sqrt_log	= u[ 0 ] * u[ 0 ] + u[ 1 ] * u[ 1 ];					/* Only points in the unit disk are accepted */
		} while( sqrt_log >= 1.0 );												/* The rejection rate is about 1 - \pi / 4 ~ 21%. */
		sqrt_log = sqrtl( -2 * logl( sqrt_log ) / sqrt_log );
		normals[ 0 ] = sqrt_log * u[ 0 ]; normals[ 1 ] = sqrt_log * u[ 1 ];
	}
	normals_ind = (normals_ind + 1) & 0x01; /* normals_ind \in {0,1}, so normals_ind := (normals_ind + 1) % 2 */
	if( !(probability( spiking_neuron ) > dsfmt_genrand_close_open( double_rngs ) * max[ spiking_neuron ]) ) {
		/* The spike is rejected */
		++nb_rejected;
		goto NEXT_SPIKE;
	}
	/* From now on the jump has been accepted,  */
	++ nb_accepted;
	if( nb_itr_per_cent >= 1 && nb_accepted % nb_itr_per_cent == 0 ) {
		fprintf( stdout, "Passed %d%% of accepted spikes.\n", nb_accepted / nb_itr_per_cent );
	}
	y_last[ spiking_neuron ] = reset_value[ spiking_neuron ];

	spiking_times[ spiking_times_array_index ] = time_int.lower_bound + delta_t;
	++spiking_times_array_index;
	if( spiking_times_array_index > size ) {
		fwrite( spiking_times, sizeof( long double ), nb_neurons, f_results );
		spiking_times_array_index = 0;
	}
	switch( conn_e ) {
		case RANDOM:			/* The indices of the postsynaptic neurons are stored in the interaction_graph array */
								for( i = 0; i < nb_couplings[ spiking_neuron ]; ++i ) {
									y_last[ interaction_graph[ spiking_neuron ][ i ] ] += interaction( spiking_neuron, i );
								}
			break;
		case FULL:			/* All neurons are postsynaptic neurons */
								for( i = 0; i < nb_neurons; ++i ) {
									y_last[ i ] += interaction( spiking_neuron, i );
								}
			break;
		case COMPLETE:			/* All neurons are postsynaptic neurons */
								for( i = 0; i < nb_neurons; ++i ) {
									if( i != spiking_neuron ){
										y_last[ i ] += interaction( spiking_neuron, i );
									}
								}
			break;
		case INDEPENDENT:		/* There are no connections between neurons, nothing to be done in this case. */
			break;
		case RECONSTRUCTION:	/* Interaction graph was too big, must be reconstructed on the fly */
								save_rng = *int_rngs;
								*int_rngs = reconstruction_graph[ spiking_neuron ];
								/* At this point it is assumed the array of indices 'indices' is sorted */
								max_ind = size;
								for( i = 0; i < nb_couplings[ spiking_neuron ]; ++i ) {
									j = get_int( max_ind ); /* A random non already chosen index in the array of indices (meta) */
									y_last[ indices[ j ] ] += interaction( spiking_neuron, j );

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
			break;
	}

	time_int.upper_bound = time_int.lower_bound + delta_t;
	goto NEXT_STEP;
END:
	if( spiking_times_array_index > 0 ) {
		fwrite( spiking_times, sizeof( long double ), spiking_times_array_index, f_results );
	}
	fclose( f_results );
}

inline
long double interaction( const unsigned int presynaptic, const unsigned int postsynaptic ) {
	(void)presynaptic; (void)postsynaptic; /* The parameters are not used currently */
	return inter;
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
	return (prob > 0)? slope * powl( prob, exponent ): 0.0;
}

inline
long double max_prob( const unsigned int neuron, const time_interval_t* const time_int ) {
	long double max_prob;
	if( fabsl( lambda[ neuron ] ) < EPSILON ) {
		max_prob = y_last[ neuron ] + var[ neuron ]
				 - threshold[ neuron ];
	} else {
		max_prob = a[ neuron ]
					+ expl( -lambda[ neuron ] * (time_int->upper_bound - t_last[ neuron ]) ) * (y_last[ neuron ] - a[ neuron ])
					+ var[ neuron ] * sqrtl(1 - expl( -2 * lambda[neuron] * (time_int->upper_bound - t_last[neuron]) ))
				 - threshold[ neuron ];
	}
	return (max_prob > 0)? slope * powl( max_prob, exponent ): 0.0;
}

