#include <stdio.h>
#include <stdlib.h>

#define N 1200

int main( int argc, char *argv[] )
{
    FILE *f;
    double **data;
    size_t i, j;

    if ( argc < 2 ) {
        printf( "usage: %s <data>\n", argv[0] );
        exit( 0 );
    }
    data = (double **) malloc( (N+1) * sizeof(double *) );
    for ( i = 0; i <= N; i++ ) {
        data[i] = (double *) malloc( (N+1) * sizeof(double) );
    }
    f = fopen( argv[1], "r" );
    for ( j = 0; j <= N; j++ ) {
        for ( i = 0; i <= N; i++ ) {
            fscanf( f, "%lf", &data[i][j] );
        }
    }
    fclose( f );
    f = fopen( argv[1], "w" );
    for ( j = 0; j <= N; j++ ) {
        for ( i = 0; i <= N; i++ ) {
            fprintf( f, "%+.16lf ", data[j][i] );
        }
        fprintf( f, "\n" );
    }
    for ( i = 0; i <= N; i++ ) {
        free( data[i] );
    }
    free( data );
    return 0;
}