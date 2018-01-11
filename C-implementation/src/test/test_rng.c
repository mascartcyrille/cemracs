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
#include "../lib/SFMT-src-1.5.1/SFMT.h"		/* Simple precision rng */
#include "../lib/dSFMT-src-2.2.3/dSFMT.h"	/* Simple precision rng */

#define	NB_RNGS	5

unsigned int	nb_neurons,
				ind, max_ind,
				size,
				i, j,  k,
				it,
				* nb_couplings,
				* indices,
				* tmp_neurons,
				** interaction_graph;
double conn_prob;
sfmt_t	*	int_rngs,
		*	tmp_int_rngs;
dsfmt_t * double_rngs;
sfmt_t * reconstruction_graph;
sfmt_t * reconstruction_graph2;
FILE * file;
char ** file_names;

unsigned int get_int2( const unsigned int limit, sfmt_t * int_rng ) {
	unsigned int ret_val, u, divide;

	divide = UINT_MAX / (limit + 1);	/* {0, ..., limit}, so limit+1 possibilities */
	do {
		u = sfmt_genrand_uint64( int_rng );
		ret_val = u / divide;
	} while( ret_val > limit );
	return ret_val;
}

unsigned int get_int( const unsigned int limit ) {
	unsigned int ret_val, u, divide;

	divide = UINT_MAX / (limit + 1);	/* {0, ..., limit}, so limit+1 possibilities */
	do {
		u = sfmt_genrand_uint64( int_rngs );
		ret_val = u / divide;
	} while( ret_val > limit );
	return ret_val;
}

int main( int argc, char ** argv ) {
/********************************************************************************/
/*									Common part									*/
/********************************************************************************/
	clock_t s11 = 0, s12 = 0, s13 = 0, s14 = 0,
			s21 = 0, s22 = 0, s23 = 0, s24 = 0,
			s31 = 0, s32 = 0, s33 = 0, s34 = 0,
			s41 = 0, s42 = 0, s43 = 0, s44 = 0,
			tmp_s;
	nb_neurons = (unsigned int) strtoul( argv[ 1 ], NULL, 10 );
	conn_prob = (double) strtod( argv[ 2 ], NULL );
	printf( "Neurons: %u\nConnection: %f\n", nb_neurons, conn_prob );
	size = nb_neurons - 1;
	(void) argc; (void) argv;
	
	tmp_neurons = malloc( sizeof( unsigned int ) * nb_neurons );
	nb_couplings = malloc( sizeof( unsigned int ) * nb_neurons );
	indices = malloc( sizeof( unsigned int ) * nb_neurons );
	int_rngs = malloc( sizeof( sfmt_t ) );
	tmp_int_rngs = malloc( sizeof( sfmt_t ) );
	double_rngs = malloc( sizeof( dsfmt_t ) );

	sfmt_init_gen_rand( int_rngs, 12345 );
	dsfmt_init_gen_rand( double_rngs, 12345 );

	for( i = 0; i < nb_neurons; ++i ) {
		/* Draw a random integer between 0 and nb_neurons included */
		nb_couplings[ i ] = 0;
		tmp_neurons[ i ] = 0;
		for( j = 0; j < nb_neurons; ++j ) {
			if( dsfmt_genrand_close_open( double_rngs ) < conn_prob ) {
				++(nb_couplings[ i ]);
			}
		}
	}

	for( i = 0; i < nb_neurons; ++i ) {
		indices[ i ] = i;
	}

// /********************************************************************************/
// /*									In memory part								*/
// /********************************************************************************/
// 	if( nb_neurons < 50000 ) {
// 		printf( "In memory part\n" );
// 		s11 = clock();
// 		interaction_graph = malloc( sizeof( unsigned int * ) * nb_neurons );
// 		for( i = 0; i < nb_neurons; ++i ) {
// 			max_ind = size;
// 			interaction_graph[ i ] = malloc( sizeof( unsigned int ) * nb_couplings[ i ] );
// 			for ( j = 0; j < nb_couplings[ i ]; ++j ) {
// 				ind = get_int( max_ind );
// 				interaction_graph[ i ][ j ] = indices[ ind ];
// 				indices[ ind ] = indices[ max_ind ];
// 				indices[ max_ind ] = interaction_graph[ i ][ j ];
// 				--max_ind;
// 			}
// 		}
// 		for( i = 0; i < nb_neurons; ++i ) {
// 			indices[ i ] = i;
// 		}
// 		s12 = clock();
// 		tmp_s = s12;
		
// 		for( it = 0; it < nb_neurons; ++it ){
// 			// Little change of the rng state
// 			i = get_int( size );
// 			// Recovering an ancient rng state
// 			s13 += clock() - tmp_s;
// 			max_ind = size;
// 			// #pragma omp parallel for
// 			for ( j = 0; j < nb_couplings[ i ]; ++j ) {
// 				tmp_neurons[ j ] = interaction_graph[ i ][ j ];
// 			}
// 			s14 = clock();
// 			tmp_s = s14;
// 		}
// 		s13 += clock() - tmp_s;
// 		for( i = 0; i < nb_neurons; ++i ) {
// 			free( interaction_graph[ i ] );
// 		}
// 		free( interaction_graph );
// 		printf( "-----\n" );
// 	}

/********************************************************************************/
/*								Reconstruction part								*/
/********************************************************************************/
	// Initialisation
	printf( "Reconstruction part\n" );
	s21 = clock();
	reconstruction_graph = malloc( sizeof( sfmt_t ) * nb_neurons );
	for( i = 0; i < nb_neurons; ++i ) {
		reconstruction_graph[ i ] = *int_rngs;
		max_ind = size;
		for ( j = 0; j < nb_couplings[ i ]; ++j ) {
			get_int( max_ind );
			--max_ind;
		}
	}
	// printf( "\n" );
	s22 = clock();
	tmp_s = s22;
	
	for( it = 0; it < nb_neurons; ++it ){
		// Little change of the rng state
		i = get_int( size );
		*tmp_int_rngs = *int_rngs;

		// Recovering an ancient rng state
		s23 += clock() - tmp_s;
		*int_rngs = reconstruction_graph[ i ];
		max_ind = size;
		for ( j = 0; j < nb_couplings[ i ]; ++j ) {
			k = get_int( max_ind );

			tmp_neurons[ j ] = i;
			ind = indices[ k ];
			indices[ k ] = indices[ max_ind ];
			indices[ max_ind ] = ind;

			--max_ind;
		}
		for( i = 0; i < nb_neurons; ++i ) {
			indices[ i ] = i;
		}
		*int_rngs = *tmp_int_rngs;
		s24 = clock();
		tmp_s = s24;
	}
	s23 += clock() - tmp_s;

	free( reconstruction_graph );
	printf( "-----\n" );

/********************************************************************************/
/*								Reconstruction part								*/
/********************************************************************************/
	// Initialisation
	free( int_rngs );
	int_rngs = malloc( sizeof( sfmt_t ) * nb_neurons );
	for( i = 0; i < nb_neurons; ++i ) {
		sfmt_init_gen_rand( &(int_rngs[ i ]), 1 + i );
	}
	printf( "Reconstruction part\n" );
	s31 = clock();
	reconstruction_graph2 = malloc( sizeof( sfmt_t ) * nb_neurons );
	// #pragma omp parallel for shared( size, reconstruction_graph2, int_rngs, nb_couplings ) private( max_ind, j ) // schedule( dynamic, 1000 )
	for( i = 0; i < nb_neurons; ++i ) {
		reconstruction_graph2[ i ] = int_rngs[ i ];
		max_ind = size;
		for ( j = 0; j < nb_couplings[ i ]; ++j ) {
			get_int2( max_ind, &(int_rngs[ i ]) );
			--max_ind;
		}
	}
	s32 = clock();
	tmp_s = s32;
	for( it = 0; it < nb_neurons; ++it ){
		// Little change of the rng state
		i = get_int2( size, int_rngs );
		// Recovering an ancient rng state
		s33 += clock() - tmp_s;
		*tmp_int_rngs = reconstruction_graph2[ i ];
		max_ind	= size;
		for( j = 0; j < nb_couplings[ i ]; ++j ) {
			k = get_int2( max_ind, tmp_int_rngs );

			tmp_neurons[ j ] = i;
			ind = indices[ k ];
			indices[ k ] = indices[ max_ind ];
			indices[ max_ind ] = ind;

			--max_ind;
		}
		for( i = 0; i < nb_neurons; ++i ) {
			indices[ i ] = i;
		}
		s34 = clock();
		tmp_s = s34;
	}
	s33 += clock() - tmp_s;

	free( reconstruction_graph2 );
	printf( "-----\n" );

// /********************************************************************************/
// /*							File part: hard drive								*/
// /********************************************************************************/
// 	// Initialisation
// 	printf( "File part: hard drive\n" );
// 	s31 = clock();
// 	file_names = malloc( sizeof( char * ) * nb_neurons );
// 	for( i = 0; i < nb_neurons; ++i ) {
// 		j = 0;
// 		file_names[ i ] = malloc( sizeof( char ) * 19 );
// 		file_names[ i ][ j++ ] = 't';
// 		file_names[ i ][ j++ ] = 'm';
// 		file_names[ i ][ j++ ] = 'p';
// 		file_names[ i ][ j++ ] = '/';
// 		file_names[ i ][ j++ ] = 't';
// 		file_names[ i ][ j++ ] = 'm';
// 		file_names[ i ][ j++ ] = 'p';
// 		file_names[ i ][ j++ ] = '_';
// 		sprintf( &(file_names[ i ][ j++ ]), "%010u", i );
// 	}
// 	for( i = 0; i < nb_neurons; ++i ) {
// 		max_ind = size;
// 		file = fopen( file_names[ i ], "wb" );
// 		if( file == NULL ) {
// 			printf( "Error while opening file %s\n", file_names[ i ] );
// 			exit( -1 );
// 		}
// 		for ( j = 0; j < nb_couplings[ i ]; ++j ) {
// 			k = get_int( max_ind );
// 			ind = indices[ k ];
// 			tmp_neurons[ j ] = ind;
// 			indices[ k ] = indices[ max_ind ];
// 			indices[ max_ind ] = ind;

// 			--max_ind;
// 		}
// 		fwrite( tmp_neurons, sizeof( unsigned int ), nb_neurons, file );
// 		fclose( file );
// 		for( j = 0; j < nb_neurons; ++j ) {
// 			indices[ j ] = j;
// 		}
// 	}
// 	s32 = clock();
// 	tmp_s = s32;

// 	for( it = 0; it < nb_neurons; ++it ) {
// 		// Little change of the rng state
// 		i = get_int( size );
// 		*tmp_int_rngs = *int_rngs;
	
// 		// Recovering an ancient rng state
// 		s33 += clock() - tmp_s;
// 		max_ind = size;
// 		file = fopen( file_names[ i ], "r" );
// 		k = fread( tmp_neurons, sizeof( unsigned int ), nb_neurons, file );
// 		for( j = 0; j < nb_couplings[ i ]; ++j ) {
// 			tmp_neurons[ j ] = i;
// 		}
// 		fclose( file );
// 		s34 = clock();
// 		tmp_s = s34;
// 	}
// 	s33 += clock() - tmp_s;

// 	for( i = 0; i < nb_neurons; ++i ) {
// 		free( file_names[ i ] );
// 	}
// 	printf( "-----\n" );

// ******************************************************************************
// /*							File part: soft drive								*/
// /********************************************************************************/
// 	// Initialisation
// 	printf( "File part: soft drive\n" );
// 	s41 = clock();
// 	file_names = malloc( sizeof( char* ) * nb_neurons );
// 	for( i = 0; i < nb_neurons; ++i ) {
// 		j = 0;
// 		file_names[ i ] = malloc( sizeof( char ) * 24 );
// 		file_names[ i ][ j++ ] = '/';
// 		file_names[ i ][ j++ ] = 'd';
// 		file_names[ i ][ j++ ] = 'e';
// 		file_names[ i ][ j++ ] = 'v';
// 		file_names[ i ][ j++ ] = '/';
// 		file_names[ i ][ j++ ] = 's';
// 		file_names[ i ][ j++ ] = 'h';
// 		file_names[ i ][ j++ ] = 'm';
// 		file_names[ i ][ j++ ] = '/';
// 		file_names[ i ][ j++ ] = 't';
// 		file_names[ i ][ j++ ] = 'm';
// 		file_names[ i ][ j++ ] = 'p';
// 		file_names[ i ][ j++ ] = '_';
// 		sprintf( &(file_names[ i ][ j++ ]), "%010u", i );
// 	}
// 	for( i = 0; i < nb_neurons; ++i ) {
// 		max_ind = size;
// 		file = fopen( file_names[ i ], "wb" );
// 		if( file == NULL ) {
// 			printf( "Error while opening file %s\n", file_names[ i ] );
// 		}
// 		for ( j = 0; j < nb_couplings[ i ]; ++j ) {
// 			k = get_int( max_ind );
// 			ind = indices[ k ];
// 			tmp_neurons[ j ] = ind;
// 			indices[ k ] = indices[ max_ind ];
// 			indices[ max_ind ] = ind;

// 			--max_ind;
// 		}
// 		fwrite( tmp_neurons, sizeof( unsigned int ), nb_neurons, file );
// 		fclose( file );
// 		for( j = 0; j < nb_neurons; ++j ) {
// 			indices[ j ] = j;
// 		}
// 	}
// 	s42 = clock();
// 	tmp_s = s42;

// 	for( it = 0; it < nb_neurons; ++it ) {
// 		// Little change of the rng state
// 		i = get_int( size );
// 		*tmp_int_rngs = *int_rngs;
	
// 		// Recovering an ancient rng state
// 		s43 += clock() - tmp_s;
// 		max_ind = size;
// 		file = fopen( file_names[ i ], "r" );
// 		k = fread( tmp_neurons, sizeof( unsigned int ), nb_neurons, file );
// 		for( j = 0; j < nb_couplings[ i ]; ++j ) {
// 			tmp_neurons[ j ] = i;
// 		}
// 		fclose( file );
// 		s44 = clock();
// 		tmp_s = s44;
// 	}
// 	s43 += clock() - tmp_s;
// 	printf( "-----\n" );

/********************************************************************************/
/*									Freeing memory								*/
/********************************************************************************/
	if( file_names != NULL )	
		for( i = 0; i < nb_neurons; ++i ) {
			free( file_names[ i ] );
		}
	
	free( file_names );
	free( nb_couplings );
	free( indices );
	free( tmp_neurons );
	free( int_rngs );
	free( tmp_int_rngs );
	free( double_rngs );
	
/********************************************************************************/
/*								Printing speed result							*/
/********************************************************************************/
	(void) s11; (void) s12; (void) s13; (void) s14;
	(void) s21; (void) s22; (void) s23; (void) s24;
	(void) s31; (void) s32; (void) s33; (void) s34;
	(void) s41; (void) s42; (void) s43; (void) s44;
	printf( "Init rec: %f\t\tTotal: %f\n", (double) (s12 - s11) / CLOCKS_PER_SEC, (double) (s12 - s11) / CLOCKS_PER_SEC + (double) s13 / CLOCKS_PER_SEC );
	printf( "Use reco: %f\n", (double) s13 / CLOCKS_PER_SEC );
	printf( "Init rec: %f\t\tTotal: %f\n", (double) (s22 - s21) / CLOCKS_PER_SEC, (double) (s22 - s21) / CLOCKS_PER_SEC + (double) s23 / CLOCKS_PER_SEC );
	printf( "Use reco: %f\n", (double) s23 / CLOCKS_PER_SEC );
	printf( "Init rec: %f\t\tTotal: %f\n", (double) (s32 - s31) / CLOCKS_PER_SEC, (double) (s32 - s31) / CLOCKS_PER_SEC + (double) s33 / CLOCKS_PER_SEC );
	printf( "Use file: %f\n", (double) s33 / CLOCKS_PER_SEC );
	printf( "Init rec: %f\t\tTotal: %f\n", (double) (s42 - s41) / CLOCKS_PER_SEC, (double) (s42 - s41) / CLOCKS_PER_SEC + (double) s43 / CLOCKS_PER_SEC );
	printf( "Use file: %f\n", (double) s43 / CLOCKS_PER_SEC );
	return 0;
}
