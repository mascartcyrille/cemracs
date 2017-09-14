#include <stdio.h>
#include <stdlib.h>

FILE *source;
FILE *destination;

int main( int argc, char const **argv ) {
	long double buffer;
	unsigned char f, i;
	char * file_name;
	
	file_name = malloc( sizeof( char ) * 100 );
	
	for( f = 1; f < argc; ++f ) {
		source = fopen( argv[ f ], "rb" );
		if (source) {
			i = 0;
			while( argv[ f ][ i ] != '\0' ) {
				file_name[ i ] = argv[ f ][ i ];
				++i;
			}
			file_name[ i ] = '\0';
			file_name[ i - 1 ] = 't';
			file_name[ i - 2 ] = 'x';
			file_name[ i - 3 ] = 't';
			destination = fopen( file_name, "w" );
			if(destination) {
				fread( &buffer, 1, sizeof( long double ), source );
				fprintf( destination, "%Lf", buffer );
				 while (!feof(source)) {
					fread( &buffer, 1, sizeof( long double ), source );
					fprintf( destination, ",\n%Lf", buffer );
				}
			} else {
				fprintf( stdout, "fail\n" );
			}
		} else {
			fprintf( stdout, "fail\n" );
		}

		fclose(source);
		fclose(destination);
	}

    return 0;
}
