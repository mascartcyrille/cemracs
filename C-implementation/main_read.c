#include <stdio.h>
#include <stdlib.h>

FILE *source;
FILE *destination;

int main( int argc, char const **argv ) {
	long double buffer;

	(void)argc;
	(void)argv;
	source = fopen("./result.bin", "rb");

	if (source) {
		destination = fopen("./result.txt", "w");
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

    return 0;
}