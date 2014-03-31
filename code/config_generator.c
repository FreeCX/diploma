#include <stdio.h>
#include <string.h>

int d[2] = { 100, 200 };
double alpha2[3] = { -0.0625, 0.1, 0.5 };
double eta[3] = { 0.0, 0.3, 1.41 };
double e[3] = { 0.7, 1.3, 1.41 };
const char fmt[] = 
	"int Nx = 100\n"
	"int Ny = 100\n"
	"int d = %d\n"
	"double alpha1 = -1.0\n"
	"double alpha2 = %.4f\n"
	"double beta1 = 1.0\n"
	"double beta2 = 0.25\n"
	"double epsilon = 0.001\n"
	"double etta = %.4f\n"
	"double ax = 0.025\n"
	"double ay = 0.025\n"
	"double e = %.4f\n";

int main( void )
{
	FILE *f;
	int counter = 0;
	int i, j, k, l;
	char buffer[512];
	char f_name[32];

	for ( i = 0; i < 2; i++ ) {
        for ( l = 0; l < 3; l++ ) {
		    for ( j = 0; j < 3; j++ ) {
			    for ( k = 0; k < 3; k++ ) {
				    sprintf( buffer, fmt, d[i], alpha2[l], eta[j], e[k] );
				    sprintf( f_name, "input_%02d.txt", counter );
				    counter++;
				    f = fopen( f_name, "w" );
				    fwrite( buffer, sizeof(char), strlen( buffer ), f );
				    fclose( f );
				    printf( "%s\n", f_name );
			    }
		    }
        }
	}
	return 0;
}
