#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <inttypes.h>
#include <limits.h>

#include "./SFMT-src-1.5.1/SFMT.h"				/* Simple precision rng */

sfmt_t	*	int_rngs;				/* Array [1..nb_threads+1]: one per thread plus one for the main thread */

void check_alloc( const void* const var, const char * const var_name ) {
	if( var == NULL ) fprintf( stderr, "Cannot allocate variable %s\n", var_name );
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
	unsigned int	*	nb_couplings;			/* Array [1..N]: number of postsynaptic neurons */
	sfmt_t			*	reconstruction_graph;	/* Array [1..N]: the state of the rng for reconstruction */
	unsigned int	*	memory;
	unsigned int		i, j,
						max_ind,
						nb_neurons;

	(void)argc;
	(void)argv;
	nb_neurons = 10;

	check_alloc( nb_couplings = malloc( sizeof( unsigned int ) * nb_neurons ), "nb_couplings" );
	for( i = 0; i < nb_neurons; ++i ) {
		nb_couplings[ i ] = get_int( nb_neurons );
	}
	check_alloc( reconstruction_graph = malloc( sizeof( sfmt_t ) * nb_neurons ), "reconstruction_graph" );
	check_alloc( memory = malloc( sizeof( unsigned int ) * nb_neurons ), "memory" );
	for( i = 0; i < nb_neurons; ++i ) {
		reconstruction_graph[ i ] = *int_rngs;
		max_ind = nb_neurons - 1;
		for ( j = 0; j < nb_couplings[ i ]; ++j ) {
			memory[ j ] = get_int( max_ind );
			--max_ind;
		}
	}

	for( i = 0; i < nb_neurons; ++i ) {
		*int_rngs = reconstruction_graph[ i ];
		max_ind = nb_neurons - 1;
		for( j = 0; j < nb_couplings[ i ]; ++j ) {
			if( memory[ j ] != get_int( max_ind ) ) {
				fprintf( stdout, "Not equal value found at index %u,%u!\n", i, j );
			}
			--max_ind;
		}
	}
	return 0;
}