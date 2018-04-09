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
				nb_threads,
				min, max, tid,
				ind, max_ind,
				size,
				i, j,  k,
				it,
				* nb_couplings,
				** indices,
				* tmp_neurons,
				** interaction_graph;
double conn_prob;
sfmt_t	*	int_rngs,
		*	tmp_int_rngs;
dsfmt_t * double_rngs;
sfmt_t ** reconstruction_graph;
FILE * file;
char ** file_names;

unsigned int get_int( const unsigned int limit, sfmt_t * int_rng ) {
	unsigned int ret_val, u, divide;

	divide = UINT_MAX / (limit + 1);	/* {0, ..., limit}, so limit+1 possibilities */
	do {
		u = sfmt_genrand_uint64( int_rng );
		ret_val = u / divide;
	} while( ret_val > limit );
	return ret_val;
}

int main( int argc, char ** argv ) {
/********************************************************************************/
/*									Common part									*/
/********************************************************************************/
	clock_t s11 = clock(), s12 = clock(), s13 = clock(), s14 = clock(),
			s21 = clock(), s22 = clock(), s23 = clock(), s24 = clock(),
			s31 = clock(), s32 = clock(), s33 = clock(), s34 = clock(),
			s41 = clock(), s42 = clock(), s43 = clock(), s44 = clock(),
			tmp_s;
	nb_neurons = (unsigned int) strtoul( argv[ 1 ], NULL, 10 );
	conn_prob = (double) strtod( argv[ 2 ], NULL );
	if( argc >= 4 ) {	
		nb_threads = 1 + (unsigned int) strtoul( argv[ 3 ], NULL, 10 );
	} else {
		nb_threads = 1;
	}
	printf( "Neurons: %u\nConnection: %f\nNumber of threads: %u\n", nb_neurons, conn_prob, nb_threads );
	size = nb_neurons - 1;
	
	tmp_neurons = malloc( sizeof( unsigned int ) * nb_neurons );
	nb_couplings = malloc( sizeof( unsigned int ) * nb_neurons );
	indices = malloc( sizeof( unsigned int * ) * nb_threads );
	int_rngs = malloc( sizeof( sfmt_t ) * nb_threads );
	tmp_int_rngs = malloc( sizeof( sfmt_t ) * nb_threads );
	double_rngs = malloc( sizeof( dsfmt_t ) * nb_threads );

	for( i = 0; i < nb_threads; ++i ) {
		sfmt_init_gen_rand( &(int_rngs[ i ]), 12345 + i );
		dsfmt_init_gen_rand( &(double_rngs[ i ]), 12345 + i );
		indices[ i ] = malloc( sizeof( unsigned int ) * nb_neurons );
		#pragma omp parallel for shared( i, indices, nb_neurons )
		for( j = 0; j < nb_neurons; ++j ) {
			indices[ i ][ j ] = j;
		}
	}
	#pragma omp parallel shared( nb_couplings, double_rngs, conn_prob, tmp_neurons ) private( tid, j ) num_threads( nb_threads )
	{
		tid = omp_get_thread_num();
		#pragma omp for
		for( i = 0; i < nb_neurons; ++i ) {
			/* Draw a random integer between 0 and nb_neurons included */
			nb_couplings[ i ] = 0;
			tmp_neurons[ i ] = 0;
			for( j = 0; j < nb_neurons; ++j ) {
				if( dsfmt_genrand_close_open( &(double_rngs[ tid ]) ) < conn_prob ) {
					++(nb_couplings[ i ]);
				}
			}
		}
	}
	// for( i = 0; i < nb_threads; ++i ) {
	// 	for( j = 0; j < nb_neurons; ++j ) {
	// 		printf( "%u ", indices[ i ][ j ] );
	// 	}
	// 	printf( "\n" );
	// }

	// for( i = 0; i < nb_neurons; ++i ) {
	// 	printf( "Neuron %u has %u children\n", i, nb_couplings[ i ] );
	// }

// /********************************************************************************/
// /*									In memory part								*/
// /********************************************************************************/
// 	if( nb_neurons < 50000 ) {
// 		printf( "In memory part\n" );
// 		s11 = clock();
// 		interaction_graph = malloc( sizeof( unsigned int * ) * nb_neurons );
// 		for( i = 0; i < nb_neurons; ++i ) {
// 			interaction_graph[ i ] = malloc( sizeof( unsigned int ) * nb_couplings[ i ] );
// 			#pragma omp parallel shared( nb_couplings, i, interaction_graph, indices, size ) private( j, max_ind, tid, ind ) num_threads( nb_threads )
// 			{
// 				max_ind = size;
// 				tid = omp_get_thread_num();
// 				#pragma omp for
// 				for ( j = 0; j < nb_couplings[ i ]; ++j ) {
// 					ind = get_int( max_ind, &(int_rngs[ tid ]) );
// 					interaction_graph[ i ][ j ] = indices[ tid ][ ind ];
// 					indices[ tid ][ ind ] = indices[ tid ][ max_ind ];
// 					indices[ tid ][ max_ind ] = interaction_graph[ i ][ j ];
// 					--max_ind;
// 				}
// 			}
// 		}
// 		#pragma omp parallel shared( indices, nb_neurons ) num_threads( nb_threads )
// 		{
// 			tid = omp_get_thread_num();
// 			#pragma omp for
// 			for( i = 0; i < nb_neurons; ++i ) {
// 				indices[ tid ][ i ] = i;
// 			}
// 		}
// 		s12 = clock();
// 		tmp_s = s12;
		
// 		for( it = 0; it < nb_neurons; ++it ){
// 			// Little change of the rng state
// 			i = get_int( size, int_rngs );
// 			// Recovering an ancient rng state
// 			s13 += clock() - tmp_s;
// 			max_ind = size;
// 			#pragma omp parallel for
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

// /********************************************************************************/
// /*								Reconstruction part								*/
// /********************************************************************************/
// 	// Initialisation
// 	printf( "Reconstruction part\n" );
// 	s21 = clock();
// 	reconstruction_graph = malloc( sizeof( sfmt_t * ) * nb_neurons );
// 	for( i = 0; i < nb_neurons; ++i ) {
// 		reconstruction_graph[ i ] = malloc( sizeof( sfmt_t ) * nb_threads );
// 	}
	
// 	for( i = 0; i < nb_neurons; ++i ) {
// 		#pragma omp parallel shared( reconstruction_graph, int_rngs, size, nb_couplings ) private( min, max, tid, j, max_ind ) num_threads( nb_threads )
// 		{
// 			tid = omp_get_thread_num();
// 			reconstruction_graph[ i ][ tid ] = int_rngs[ tid ];
// 			max_ind = size;
// 			min = tid * nb_couplings[ i ] / nb_threads;
// 			max	= (tid + 1) * nb_couplings[ i ] / nb_threads;
// 			#pragma omp for
// 			for ( j = min; j < max; ++j ) {
// 				get_int( max_ind, &(int_rngs[ tid ]) );
// 				--max_ind;
// 			}
// 		}
// 	}
// 	s22 = clock();
// 	tmp_s = s22;
	
// 	for( it = 0; it < nb_neurons; ++it ){
// 		// Little change of the rng state
// 		i = get_int( size, int_rngs );
// 		for( j = 0; j < nb_threads; ++j ) {
// 			tmp_int_rngs[ j ] = int_rngs[ j ];
// 		}
// 		// Recovering an ancient rng state
// 		s23 += clock() - tmp_s;
// 		#pragma omp parallel shared( reconstruction_graph, int_rngs, size, nb_couplings, indices ) private( min, max, tid, k, j, max_ind ) num_threads( nb_threads )
// 		{
// 			tid = omp_get_thread_num();
// 			min = tid * nb_couplings[ i ] / nb_threads;
// 			max	= (tid + 1) * nb_couplings[ i ] / nb_threads;
// 			int_rngs[ tid ] = reconstruction_graph[ i ][ tid ];
// 			max_ind	= size;
// 			#pragma omp for
// 			for( j = min; j < max; ++j ) {
// 				k = get_int( max_ind, &(int_rngs[ tid ]) );

// 				tmp_neurons[ j ] = i;
// 				ind = indices[ tid ][ k ];
// 				indices[ tid ][ k ] = indices[ tid ][ max_ind ];
// 				indices[ tid ][ max_ind ] = ind;

// 				--max_ind;
// 			}
// 		}
// 		#pragma omp parallel shared( indices, nb_neurons ) num_threads( nb_threads )
// 		{
// 			tid = omp_get_thread_num();
// 			#pragma omp for
// 			for( i = 0; i < nb_neurons; ++i ) {
// 				indices[ tid ][ i ] = i;
// 			}
// 		}
// 		for( i = 0; i < nb_threads; ++i ) {
// 			int_rngs[ i ] = tmp_int_rngs[ i ];
// 		}
// 		s24 = clock();
// 		tmp_s = s24;
// 	}
// 	s23 += clock() - tmp_s;

// 	for( i = 0; i < nb_neurons; ++i ) {
// 		free( reconstruction_graph[ i ] );
// 	}
// 	free( reconstruction_graph );
// 	printf( "-----\n" );

/********************************************************************************/
/*								Reconstruction part								*/
/********************************************************************************/
	// Initialisation
	free( int_rngs );
	free( tmp_int_rngs );
	int_rngs = malloc( sizeof( sfmt_t ) * nb_neurons );
	tmp_int_rngs = malloc( sizeof( sfmt_t ) );
	for( i = 0; i < nb_neurons; ++i ) {
		sfmt_init_gen_rand( &(int_rngs[ i ]), 1 + i );
	}
	printf( "Reconstruction part\n" );
	s31 = clock();
	reconstruction_graph = malloc( sizeof( sfmt_t * ) * nb_neurons );
	for( i = 0; i < nb_neurons; ++i ) {
		reconstruction_graph[ i ] = malloc( sizeof( sfmt_t ) );
	}
	// #pragma omp parallel for shared( size, reconstruction_graph, int_rngs, nb_couplings ) private( max_ind, j ) // schedule( dynamic, 1000 )
	for( i = 0; i < nb_neurons; ++i ) {
		*reconstruction_graph[ i ] = int_rngs[ i ];
		max_ind = size;
		for ( j = 0; j < nb_couplings[ i ]; ++j ) {
			get_int( max_ind, &(int_rngs[ i ]) );
			--max_ind;
		}
	}
	s32 = clock();
	tmp_s = s32;
	
	for( it = 0; it < nb_neurons; ++it ){
		// Little change of the rng state
		i = get_int( size, int_rngs );
		*tmp_int_rngs = int_rngs[ i ];
		// Recovering an ancient rng state
		s33 += clock() - tmp_s;
		*tmp_int_rngs = *reconstruction_graph[ i ];
		max_ind	= size;
		for( j = 0; j < nb_couplings[ i ]; ++j ) {
			k = get_int( max_ind, tmp_int_rngs );

			tmp_neurons[ j ] = i;
			ind = indices[ 0 ][ k ];
			indices[ 0 ][ k ] = indices[ 0 ][ max_ind ];
			indices[ 0 ][ max_ind ] = ind;

			--max_ind;
		}
		for( i = 0; i < nb_neurons; ++i ) {
			indices[ 0 ][ i ] = i;
		}
		s34 = clock();
		tmp_s = s34;
	}
	s33 += clock() - tmp_s;

	for( i = 0; i < nb_neurons; ++i ) {
		free( reconstruction_graph[ i ] );
	}
	free( reconstruction_graph );
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
// 		#pragma omp parallel for
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
// 		#pragma omp parallel for
// 		for( j = 0; j < nb_couplings[ i ]; ++j ) {
// 			tmp_neurons[ j ] = i;
// 		}
// 		fclose( file );
// 		s34 = clock();
// 		tmp_s = s34;
// 	}
// 	s33 += clock() - tmp_s;

// 	#pragma omp parallel for
// 	for( i = 0; i < nb_neurons; ++i ) {
// 		free( file_names[ i ] );
// 	}
// 	printf( "-----\n" );

// /********************************************************************************/
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
// 		#pragma omp parallel for
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
// 		#pragma omp parallel for
// 		for( j = 0; j < nb_couplings[ i ]; ++j ) {
// 			tmp_neurons[ j ] = i;
// 		}
// 		fclose( file );
// 		s44 = clock();
// 		tmp_s = s44;
// 	}
// 	s43 += clock() - tmp_s;
// 	printf( "-----\n" );
// 	for( i = 0; i < nb_neurons; ++i ) {
// 		free( file_names[ i ] );
//  }
//  free( file_names );

/********************************************************************************/
/*									Freeing memory								*/
/********************************************************************************/
	free( nb_couplings );
	for( i = 0; i < nb_threads; ++i ) {
		free( indices[ i ] );
	}
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
	(void) tmp_s;
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
