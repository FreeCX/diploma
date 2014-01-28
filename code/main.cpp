#include <cstdio>
#include <cmath>
#include <GL/glut.h>
#include "alglib/specialfunctions.h"

const float q0 = 1.0, q1 = 1.1, q2 = 0.9;
const float m0 = 1.0, m1 = 2.0, m2 = sqrt( 2.0 );

const float range = 2.00f;
const float delta = 0.05f;
const unsigned long N = pow( (unsigned long)( 2*range / delta ) + 1, 2 );

int w_width = 800;
int w_height = 800;
int win_id = 0;
GLfloat n = 1.2f;
GLfloat aspect;

float *vertex = NULL;
float *color = NULL;

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
	vertex = new float [3*N];
	color = new float [3*N];
	for ( float x = -range; x <= range; x += delta ) {
    	for ( float y = -range; y <= range; y += delta ) {
    		vertex[3*i+0] = x;
    		vertex[3*i+1] = y;
    		vertex[3*i+2] = V( sqrt( x*x+y*y ) );
    		if ( vertex[3*i+2] < 0 ) {
    			color[3*i+0] = 0.0f;
    			color[3*i+1] = 0.0f;
    			color[3*i+2] = 1.0f;
    		} else if ( vertex[3*i+2] > 0 ) {
    			color[3*i+0] = 1.0f;
    			color[3*i+1] = 0.0f;
    			color[3*i+2] = 0.0f;
    		}
    		i++;
    	}
    }
}

void program_init( void )
{
    GLint sw = glutGet( GLUT_SCREEN_WIDTH );
    GLint sh = glutGet( GLUT_SCREEN_HEIGHT );
    glutPositionWindow( ( sw - w_width ) / 2, ( sh - w_height) / 2 );
    glClearColor( 0.0f, 0.0f, 0.0f, 1.0f );
    glPointSize( 2.0f );
    glEnable( GL_DEPTH_TEST );
    calculate();
}

void program_free( void )
{
	delete color;
	delete vertex;
}

void program_render( void )
{
	static float t = -2*M_PI;
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glLoadIdentity();
    glScalef( 0.5f, 0.5f, 0.5f );
    glColor3f( 1.0f, 1.0f, 1.0f );
    glRotatef( 45.0f, -1.0f, 0.0f, 0.0f );
    glRotatef( 100*t, 0.0f, 0.0f, 1.0f );
    glEnableClientState( GL_VERTEX_ARRAY );
    glEnableClientState( GL_COLOR_ARRAY );
    glVertexPointer( 3, GL_FLOAT, 0, vertex );
    glColorPointer( 3, GL_FLOAT, 0, color );
    glDrawArrays( GL_POINTS, 0, N );
    glDisableClientState( GL_COLOR_ARRAY );
    glDisableClientState( GL_VERTEX_ARRAY );
    glutSwapBuffers();
    if ( t > 2*M_PI ) {
    	t = -2*M_PI;
    } else {
    	t += 0.01f;
    }
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
    gluPerspective( 0.0f, aspect, 0.0f, 20.0f );
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
}
void program_keyboard( unsigned char key, int x, int y )
{
    static int fullscreen = 1;
    if ( key == 'q' ) {
        glutDestroyWindow( win_id );
    } else if ( key == 'f' ) {
        if ( fullscreen ) {
            glutFullScreen();
        } else {
            glutReshapeWindow( w_width, w_height );
        }
        fullscreen = !fullscreen;
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
