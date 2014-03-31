#include "parser.h"

struct u_type {
    char type[16];
    char utype;
};
typedef struct u_type u_type_t;

u_type_t t[4] = {
    { "int", T_INT },
    { "float", T_FLOAT },
    { "double", T_DOUBLE },
    { 0 }
};

int config_parser( const char *filename, int count, p_block_t *a )
{
    char *p, utype, line[1024];
    FILE *f;
    int i;

    f = fopen( filename, "r" );
    if ( f == NULL ) {
        return E_FAILED;
    }
    while ( !feof( f ) ) {
        fgets( line, 1024, f );
        p = strtok( line, " " );
        for ( i = 0; t[i].type != NULL; i++ ) {
            if ( strcmp( t[i].type, p ) == 0 ) {
                utype = t[i].utype;
                break;
            }
        }
        p = strtok( NULL, " " );
        if ( p == NULL ) {
            break;
        }
        for ( i = 0; a[i].name != NULL; i++ ) {
            if ( strcmp( a[i].name, p ) == 0 ) {
                break;
            }
        }
		if (i > count) {
			continue;
		}
        a[i].utype = (char) utype;
        p = strtok( NULL, "=" );
        switch ( a[i].utype ) {
            case T_INT:
                *a[i].ivalue = strtol( p, (char **) NULL, 10 );
                break;
            case T_FLOAT:
                *a[i].fvalue = strtof( p, (char **) NULL );
                break;
            case T_DOUBLE:
                *a[i].dvalue = strtod( p, (char **) NULL );
                break;
        }
    }
    fclose( f );
    return E_SUCCESS;
}

/* how to use block */
/*
    int param1;
    float param2;
    double param3;
    p_block_t blk[] = {
        { "param1", T_INT, &param1 },
        { "param2", T_FLOAT, 0, &param2 },
        { "param3", T_DOUBLE, 0, 0, &param3 },
        { 0 }
    };
    config_parser( "file.txt", 3, blk );
    printf( "value1 = %d\n", param1 );
    printf( "value2 = %f\n", param2 );
    printf( "value3 = %E\n", param3 );
*/