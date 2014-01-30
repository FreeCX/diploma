#include <cstdio>
#include <cmath>
#include <GL/glut.h>
#include "alglib/specialfunctions.h"

float q0 = 1.0, q1 = 1.1, q2 = 0.9;
float m0 = 1.0, m1 = 2.0, m2 = sqrt( 2.0 );

float adder = 0.1f;
int state = 0;

const float range = 2.00f;
const float delta = 0.04f;
const unsigned long N = pow( (unsigned long)( 2*range / delta ) + 1, 2 );

int w_width = 800;
int w_height = 800;
int win_id = 0;
GLfloat n = 1.2f;
GLfloat aspect;

float *vertex = NULL;
float *color = NULL;

void renderFont( float x, float y, void *font, const char *fmt, ... ) 
{
    char text[256], *c;
    va_list ap;
    float x1 = x;
    
    va_start( ap, fmt );
    vsprintf( text, fmt, ap );
    va_end( ap );
    for ( c = text; *c != '\0'; c++ ) {
        glRasterPos2f( x1, y );
        glutBitmapCharacter( font, *c );
        x1 += 0.07f / n;
    }
}

double V( double r )
{
    return 2.0 * M_PI * ( 
        pow( q0, 2.0 ) * alglib::besselk0( m0 * r ) - 
        pow( q1, 2.0 ) * alglib::besselk0( m1 * r ) - 
        pow( q2, 2.0 ) * alglib::besselk0( m2 * r ) 
    );
}

void calculate( void )
{
    int i = 0;
    float min, max, mid;
    if ( vertex == NULL ) {
        vertex = new float [3*N];
        color = new float [3*N];
    }
    min = max = V( sqrt( 2*range*range ) );
    for ( float x = -range; x <= range; x += delta ) {
        for ( float y = -range; y <= range; y += delta ) {
            vertex[3*i+0] = x;
            vertex[3*i+1] = y;
            vertex[3*i+2] = V( sqrt( x*x+y*y ) );
            if ( min > vertex[3*i+2] ) {
                min = vertex[3*i+2];
            }
            if ( max < vertex[3*i+2] ) {
                max = vertex[3*i+2];
            }
            i++;
        }
    }
    mid = 0.5 * ( max + min );
    for ( int i = 0; i < N; i++ ) {
        if ( vertex[3*i+2] > mid ) {
            color[3*i+0] = max / vertex[3*i+2];
            color[3*i+2] = 0.0f;
        } else {
            color[3*i+0] = 0.0f;
            color[3*i+2] = min / vertex[3*i+2];
        }
        color[3*i+1] = mid / vertex[3*i+2];
        vertex[3*i+2] = 0.0f;
    }
}

void program_init( void )
{
    GLint sw = glutGet( GLUT_SCREEN_WIDTH );
    GLint sh = glutGet( GLUT_SCREEN_HEIGHT );
    glutPositionWindow( ( sw - w_width ) / 2, ( sh - w_height) / 2 );
    glClearColor( 0.0f, 0.0f, 0.0f, 1.0f );
    glPointSize( 5.0f );
    calculate();
}

void program_free( void )
{
    delete color;
    delete vertex;
}

void program_render( void )
{
    glClear( GL_COLOR_BUFFER_BIT );
    glLoadIdentity();
    glScalef( 0.5f, 0.5f, 0.5f );
    glEnableClientState( GL_VERTEX_ARRAY );
    glEnableClientState( GL_COLOR_ARRAY );
    glVertexPointer( 3, GL_FLOAT, 0, vertex );
    glColorPointer( 3, GL_FLOAT, 0, color );
    glDrawArrays( GL_POINTS, 0, N );
    glDisableClientState( GL_COLOR_ARRAY );
    glDisableClientState( GL_VERTEX_ARRAY );
    glColor3f( 1.0f, 1.0f, 1.0f );
    renderFont( -2.3f, 2.2f, GLUT_BITMAP_8_BY_13, 
        "[%d]: q0 = %.3f, q1 = %.3f, q2 = %.3f, a = %.0E", 
        state+1, q0, q1, q2, adder );
    glutSwapBuffers();
}

void program_redraw( int value )
{
    program_render();
    glutTimerFunc( 30, program_redraw, 0 );
}

void program_resize( int width, int height )
{
    aspect = (float) width / height;
    glViewport( 0, 0, width, height );
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    if ( width <= height ) {
        glOrtho( -n, n, -n/aspect, n/aspect, n, -n );
    } else {
        glOrtho( -n * aspect, n * aspect, -n, n, n, -n );
    }
    gluPerspective( 0.0f, aspect, 0.0f, 10.0f );
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
}
void program_keyboard( unsigned char key, int x, int y )
{
    static int fullscreen = 1;
    switch ( key ) {
        case 'q':
            glutDestroyWindow( win_id );
            break;
        case 'f':
            if ( fullscreen ) {
                glutFullScreen();
            } else {
                glutReshapeWindow( w_width, w_height );
            }
            fullscreen = !fullscreen;
            break;
        case '+':
            switch ( state ) {
                case 0:
                    q0 += adder;
                    break;
                case 1:
                    q1 += adder;
                    break;
                case 2:
                    q2 += adder;
                    break;
                case 3:
                    adder *= 10.0f;
                    break;
            }
            calculate();
            break;
        case '-':
            switch ( state ) {
                case 0:
                    q0 -= adder;
                    break;
                case 1:
                    q1 -= adder;
                    break;
                case 2:
                    q2 -= adder;
                    break;
                case 3:
                    adder /= 10.0f;
                    break;
            }
            calculate();
            break;
        case '1' ... '4':
            state = key-'0'-1;
            break;
    }
}
int main( int argc, char *argv[] )
{
    glutInit( &argc, argv );
    glutInitWindowSize( w_width, w_height );
    glutInitWindowPosition( 0, 0 );
    glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE );
    win_id = glutCreateWindow( "Superconductive" );
    glutReshapeFunc( program_resize );
    glutDisplayFunc( program_render );
    glutKeyboardFunc( program_keyboard );
    glutTimerFunc( 1, program_redraw, 0 );
    program_init();
    glutMainLoop();
    program_free();
    return EXIT_SUCCESS;
}
