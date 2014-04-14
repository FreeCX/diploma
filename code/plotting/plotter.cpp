/* ------------------------------------------ */
/*   Автоматизированное построение графиков   */
/* программа: Light Plotter                   */
/* автор: dr.FreeCX                           */
/* релиз: 11/04/2014                          */
/* версия: 0.3                                */
/* ------------------------------------------ */
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>

enum {
    PLOT_3D,
    PLOT_MAP,
    PLOT_PDF,
    PLOT_PNG
};

int Nx, Ny, step;
double **data;
const char *fieldB_file = "field_B.txt";
const char *fieldF_file = "field_F%d.txt";

void create_plot( const char *input, const char *output, int type, int format )
{
    FILE *p;
    char f_format[16], f_type[16];

    p = popen( "gnuplot -p", "w" );
    fprintf( p, "set xrange[0:%d]\n", Nx );
    fprintf( p, "set yrange[0:%d]\n", Ny );
    fprintf( p, "set size square 1.0,1.0\n" );
    fprintf( p, "set palette rgbformulae 33,13,10\n" );
    printf( "  >> Plot '%s' file to '%s' .... ", input, output );
    switch ( format ) {
        case PLOT_PDF:
            fprintf( p, "set terminal pdf size 5.0in,5.0in\n" );
            strcpy( f_format, "pdf" );
            break;
        case PLOT_PNG:
            fprintf( p, "set terminal png size 800,800\n" );
            strcpy( f_format, "png" );
            break;
    }
    switch ( type ) {
        case PLOT_3D:
            fprintf( p, "set pm3d\n" );
            // fprintf( p, "set view 45,45,1,1\n" );
            strcpy( f_type, "3d" );
            break;
        case PLOT_MAP:
            fprintf( p, "set pm3d map\n" );
            strcpy( f_type, "map" );
            break;
    }
    fprintf( p, "set output '%s_%s.%s'\n", output, f_type, f_format );
    fprintf( p, "splot '%s' every %d:%d matrix "
                        "w pm3d palette title ''\n", input, step, step );
    puts( "[ Done ]" );
    pclose( p );
}

int convert_fplot( int N, const char *output_file )
{
    FILE *f0, *f1;
    std::complex<double> tmp;
    double x, y;
    char buffer[64];

    sprintf( buffer, "dataF%dr.dat", N );
    f0 = fopen( buffer, "r" );
    sprintf( buffer, "dataF%dc.dat", N );
    f1 = fopen( buffer, "r" );
    if ( f0 == NULL || f1 == NULL ) {
        printf( "  >> Can't open dataF%d{c,r} files!\n", N );
        return 0;
    }
    printf( "  >> Converting complex to abs 'data_F%d' .... ", N );
    for ( int j = 0; j <= Ny; j++ ) {
        for ( int i = 0; i <= Nx; i++ ) {
            fscanf( f0, "%lf", &x );
            fscanf( f1, "%lf", &y );
            tmp = {x, y};
            data[j][i] = pow( std::abs( tmp ), 2 );
        }
    }
    puts( "[ Done ]" );
    fclose( f0 );
    fclose( f1 );
    printf( "  >> Write data to '%s' file .... ", output_file );
    f0 = fopen( output_file, "w" );
    if ( f0 == NULL ) {
        printf( "[FAILED]\n  [-] Can't open '%s' file!\n", output_file );
        exit( 0 );
    }
    for ( int i = 0; i <= Ny; i++ ) {
        for ( int j = 0; j <= Nx; j++ ) {
            fprintf( f0, "%+.16lf ", data[j][i] );
        }
        fprintf( f0, "\n" );
    }
    fclose( f0 );
    puts( "[ Done ]" );
    return 1;
}

// int convert_plot( const char *input_file, const char *output_file )
// {
//     FILE *f;
//     int i, j, k, m;
//     double tmp;
//     char buffer[256];

//     printf( "  >> Loading '%s' file .... ", input_file );
//     f = fopen( input_file, "r" );
//     if ( f == NULL ) {
//         printf( "[FAILED]\n  [-] Can't open '%s' file!\n", input_file );
//         exit( 0 );
//     }
//     for ( int j = 0; j <= Ny; j++ ) {
//         for ( int i = 0; i <= Nx; i++ ) {
//             fscanf( f, "%lf", &data[j][i] );
//         }
//     }
//     fclose( f );
//     puts( "[ Done ]" );
//     printf( "  >> Convert '%s' file .... ", input_file );
//     f = fopen( output_file, "w" );
//     if ( f == NULL ) {
//         printf( "[FAILED]\n  [-] Can't create '%s' file!", output_file );
//         return 0;
//     }
//     for ( int i = 0; i <= Ny; i++ ) {
//         for ( int j = 0; j <= Nx; j++ ) {
//             fprintf( f, "%+.16lf ", data[j][i] );
//         }
//         fprintf( f, "\n" );
//     }
//     fclose( f );
//     puts( "[ Done ]" );
//     return 1;
// }

void clean_folder( void )
{
    const char *fmt_del = ".txt";
    struct dirent *dp;
    DIR *dir;

    dir = opendir( "." );
    if ( dir == NULL ) {
        puts( "  [-] Can't open folder ... [ Failed ]" );
        return;
    }
    while ( ( dp = readdir( dir ) ) != NULL ) {
        if ( strstr( dp->d_name, fmt_del ) != NULL ) {
            printf( "  >> Delete '%s' file .... [ %s ]\n", dp->d_name, 
                unlink( dp->d_name ) == 0 ? "Done" : "Failed" );
        }
    }
}

void main_procedure( int output_format )
{
    char input[64], output[64];

    // if ( convert_plot( "dataB.dat", fieldB_file ) ) {
    //     create_plot( fieldB_file, "field_B", PLOT_3D, output_format );
    //     create_plot( fieldB_file, "field_B", PLOT_MAP, output_format );
    // }
    create_plot( "dataB.dat", "field_B", PLOT_3D, output_format );
    create_plot( "dataB.dat", "field_B", PLOT_MAP, output_format );
    for ( int i = 1; i <= 2; i++ ) {
        sprintf( input, fieldF_file, i );
        sprintf( output, "field_F%d", i );
        convert_fplot( i, input );
        create_plot( input, output, PLOT_3D, output_format );
        create_plot( input, output, PLOT_MAP, output_format );
    }
    clean_folder();
}

int folder_scan( int output_format )
{
    struct dirent *dp;
    struct stat fst;
    DIR *dir;
    int status = 1;

    dir = opendir( "." );
    if ( dir == NULL ) {
        return !status;
    }
    while ( ( dp = readdir( dir ) ) != NULL ) {
        if ( !strcmp( dp->d_name, "." ) || !strcmp( dp->d_name, ".." ) ) {
            continue;
        }
        stat( dp->d_name, &fst );
        if ( ( fst.st_mode & S_IFMT ) == S_IFDIR ) {
            printf( ">> Enter to '%s' folder\n", dp->d_name );
            chdir( dp->d_name );
            main_procedure( output_format );
            chdir( ".." );
            printf( ">> Work in folder '%s' complete!\n", dp->d_name );
        }
    }
    closedir( dir );
    return status;
}

void init_memory( void )
{
    data = (double **) malloc( (Nx+1) * sizeof(double *) );
    for ( int i = 0; i <= Nx; i++ ) {
        data[i] = (double *) malloc( (Ny+1) * sizeof(double) );
    }
}

void clean_memory( void )
{
    for ( int i = 0; i <= Nx; i++ ) {
        free( data[i] );
    }
    free( data );
}

int main( void )
{
    Nx = Ny = 1200;
    step = 1;
    init_memory();
    folder_scan( PLOT_PNG );
    clean_memory();
    return 0;
}