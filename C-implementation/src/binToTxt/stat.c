#include <stdio.h>
#include <stdlib.h>

int main( int argc, char ** argv ) {
	long double	curr,			/*  */
				prev,			/*  */
				min,			/*  */
				max,			/*  */
				avg,			/*  */
				tmp;			/*  */
	int			nb_neurons,		/*  */
				f;				/*  */
	FILE	*	f_input;		/*  */
	
	nb_neurons = strtoul( argv[ 1 ], NULL, 10 );
	avg = 0.0;
	min = 1000;
	max = -1000;
	
	for( f = 2; f < argc; ++f ) {
		f_input = fopen( argv[ f ], "rb" );
		if( f_input ) {
			fread( &curr, 1, sizeof( long double ), f_input );
			while( !feof( f_input ) ) {
			 	prev = curr;
				fread( &curr, 1, sizeof( long double ), f_input );
				tmp = 1 / (curr - prev);
				min = ( tmp < min )? tmp: min;
				max = ( tmp > max )? tmp: max;
				avg += tmp;
				if( tmp < 0 ) {
					fprintf( stderr, "Error, end of the world, beep bip bop!\n" );
				}
			}
		} else {
			fprintf( stderr, "fail\n" );
		}
		
		fprintf( stdout, "Max: %1.16Lf\tMin: %1.16Lf\tAvg: %1.16Lf\n", max, min, avg / nb_neurons );

		fclose( f_input );
	}
	
	return 0;
}
