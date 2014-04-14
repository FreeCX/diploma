#include <complex>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <sstream>
#include <stdint.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <vector>
#include "include/alglib/src/optimization.h"
#include "cppad/cppad.hpp"
#include "HLBFGS/HLBFGS.h"
#include "parser.h"

#define USE_OPENMP 1

using namespace alglib;
using namespace std;
using CppAD::AD;

uint32_t nextTime, predTime;

// Параметры задачи
double aa = 0.0;
double alpha1 = -1;
double alpha2 = -0.0069;
double ax = 0.025;
double ay = 0.025;
double bb = 0.0;
double beta1 = 1;
double beta2 = 0.0278;
double cc = 0.0;
double deltap = 0.000001;
double e = 0.7;
double en1, en2;
double epsilon = 0.001;
double etta = 0.0111;
double F = 0.0;
double F0 = 0.0;
double F1 = 0.0;
double f1abs = 0;
double f1imag = 0;
double f1real = 0;
double F2 = 0.0;
double f2abs = 0;
double f2imag = 0;
double f2real = 0;
double H20 = 0.0;
double ksi1 = 1 / sqrt(2 * abs(alpha1));
double ksi2 = 1 / sqrt(2 * abs(alpha2));
double rax = 1 / ax;
double ray = 1 / ay;
double rdeltap = 1 / deltap;
double T0 = 0.0;
double Theta1 = 0;
double theta1, theta2, r1, r2;
double Theta2 = 0;
double u10 = sqrt( abs( alpha1 / beta1 ) );
double u20 = sqrt( abs( alpha2 / beta2 ) );
double W0 = 0.0;
double y;
int d = 1;
int iv1, jv1, iv2, jv2;
int k = 0;
int n = 0;
int Nx = 1200;
int Ny = 800;
int repeat = 1;
std::complex<double> f1 = 0, f2 = 0;
std::complex<double> I(0.0, 1.0), z1 = 0, z2 = 0;
vector< AD<double> > Y;
vector< AD<double> > X;
CppAD::ADFun<double> fun;
double **Axv, **Axy, **W, **T, **H2;
double **dAx, **dAy;
double **df1re, **df2re, **df1im, **df2im;
double *grad, *vars;
std::complex<double> *f1v, *f2v;
real_1d_array x;
// Вспомогательные переменные

void get_time( void )
{
	struct tm * ti;
	time_t raw;
	char buffer[64];

	time( &raw );
	ti = localtime( &raw );
	strftime( buffer, 64, "%d/%m/%y %H:%M:%S", ti );
	puts( buffer );
}

uint32_t get_ticks( void )
{
    struct timeval tv;
    gettimeofday( &tv, 0 );
    return ( tv.tv_sec * 1000 + tv.tv_usec / 1000 );
}

// Функция вычисления "потенциальной" энергии, связанной с узлом с координатами i, j
inline double CalculateW( int i, int j, double * vars )
{
	double f1real = vars[(Ny+1)*i+j];
	double f1imag = vars[(Ny+1)*i+j+(Nx+1)*(Ny+1)];
	double f2real = vars[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)];
	double f2imag = vars[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)];
	double Theta1 = atan2(f1imag,f1real);
	double Theta2 = atan2(f2imag,f2real);
	double f1abs = f1real*f1real+f1imag*f1imag;
	double f2abs = f2real*f2real+f2imag*f2imag;
	double W = 0.5*((2*alpha1 + beta1*f1abs)*f1abs + (2*alpha2 + beta2*f2abs)*f2abs) - etta*sqrt(f1abs*f2abs)*cos(Theta2 - Theta1);
	return W;
}

// Функция вычисления "кинетической" энергии, связанной с узлом с координатами i, j
inline double CalculateT( int i, int j, double * vars )
{
	double f1xre;
	double f2xre;
	double f1yre;
	double f2yre;
	double f1xim;
	double f2xim;
	double f1yim;
	double f2yim;

	if ( i < Nx && j < Ny && i > 0 && j > 0 ) {
		f1xre = (vars[(Ny+1)*(i+1)+j]-vars[(Ny+1)*(i-1)+j])/2*rax;
		f1xim = (vars[(Ny+1)*(i+1)+j+(Nx+1)*(Ny+1)]-vars[(Ny+1)*(i-1)+j+(Nx+1)*(Ny+1)])/2*rax;
		f2xre = (vars[(Ny+1)*(i+1)+j+2*(Nx+1)*(Ny+1)]-vars[(Ny+1)*(i-1)+j+2*(Nx+1)*(Ny+1)])/2*rax;
		f2xim = (vars[(Ny+1)*(i+1)+j+3*(Nx+1)*(Ny+1)]-vars[(Ny+1)*(i-1)+j+3*(Nx+1)*(Ny+1)])/2*rax;

		f1yre = (vars[(Ny+1)*i+j+1]-vars[(Ny+1)*i+j-1])/2*ray;
		f1yim = (vars[(Ny+1)*i+j+1+(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j-1+(Nx+1)*(Ny+1)])/2*ray;
		f2yre = (vars[(Ny+1)*i+j+1+2*(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j-1+2*(Nx+1)*(Ny+1)])/2*ray;
		f2yim = (vars[(Ny+1)*i+j+1+3*(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j-1+3*(Nx+1)*(Ny+1)])/2*ray;
	} else if ( i == 0 && j == 0 ) {
		f1xre = (vars[(Ny+1)*(i+1)+j]-vars[(Ny+1)*i+j])*rax;
		f1xim = (vars[(Ny+1)*(i+1)+j+(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j+(Nx+1)*(Ny+1)])*rax;
		f2xre = (vars[(Ny+1)*(i+1)+j+2*(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)])*rax;
		f2xim = (vars[(Ny+1)*(i+1)+j+3*(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)])*rax;

		f1yre = (vars[(Ny+1)*i+j+1]-vars[(Ny+1)*i+j])*ray;
		f1yim = (vars[(Ny+1)*i+j+1+(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j+(Nx+1)*(Ny+1)])*ray;
		f2yre = (vars[(Ny+1)*i+j+1+2*(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)])*ray;
		f2yim = (vars[(Ny+1)*i+j+1+3*(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)])*ray;
	} else if ( i == 0 && j > 0 && j < Ny ) {
		f1xre = (vars[(Ny+1)*(i+1)+j]-vars[(Ny+1)*i+j])*rax;
		f1xim = (vars[(Ny+1)*(i+1)+j+(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j+(Nx+1)*(Ny+1)])*rax;
		f2xre = (vars[(Ny+1)*(i+1)+j+2*(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)])*rax;
		f2xim = (vars[(Ny+1)*(i+1)+j+3*(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)])*rax;

		f1yre = (vars[(Ny+1)*i+j+1]-vars[(Ny+1)*i+j-1])/2*ray;
		f1yim = (vars[(Ny+1)*i+j+1+(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j-1+(Nx+1)*(Ny+1)])/2*ray;
		f2yre = (vars[(Ny+1)*i+j+1+2*(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j-1+2*(Nx+1)*(Ny+1)])/2*ray;
		f2yim = (vars[(Ny+1)*i+j+1+3*(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j-1+3*(Nx+1)*(Ny+1)])/2*ray;
	} else if ( i == 0 && j == Ny ) {
		f1xre = (vars[(Ny+1)*(i+1)+j]-vars[(Ny+1)*i+j])*rax;
		f1xim = (vars[(Ny+1)*(i+1)+j+(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j+(Nx+1)*(Ny+1)])*rax;
		f2xre = (vars[(Ny+1)*(i+1)+j+2*(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)])*rax;
		f2xim = (vars[(Ny+1)*(i+1)+j+3*(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)])*rax;

		f1yre = (vars[(Ny+1)*i+j]-vars[(Ny+1)*i+j-1])*ray;
		f1yim = (vars[(Ny+1)*i+j+(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j-1+(Nx+1)*(Ny+1)])*ray;
		f2yre = (vars[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j-1+2*(Nx+1)*(Ny+1)])*ray;
		f2yim = (vars[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j-1+3*(Nx+1)*(Ny+1)])*ray;
	} else if ( i > 0 && i < Nx && j == 0 ) {
		f1xre = (vars[(Ny+1)*(i+1)+j]-vars[(Ny+1)*(i-1)+j])/2*rax;
		f1xim = (vars[(Ny+1)*(i+1)+j+(Nx+1)*(Ny+1)]-vars[(Ny+1)*(i-1)+j+(Nx+1)*(Ny+1)])/2*rax;
		f2xre = (vars[(Ny+1)*(i+1)+j+2*(Nx+1)*(Ny+1)]-vars[(Ny+1)*(i-1)+j+2*(Nx+1)*(Ny+1)])/2*rax;
		f2xim = (vars[(Ny+1)*(i+1)+j+3*(Nx+1)*(Ny+1)]-vars[(Ny+1)*(i-1)+j+3*(Nx+1)*(Ny+1)])/2*rax;

		f1yre = (vars[(Ny+1)*i+j+1]-vars[(Ny+1)*i+j])*ray;
		f1yim = (vars[(Ny+1)*i+j+1+(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j+(Nx+1)*(Ny+1)])*ray;
		f2yre = (vars[(Ny+1)*i+j+1+2*(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)])*ray;
		f2yim = (vars[(Ny+1)*i+j+1+3*(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)])*ray;
	} else if ( i > 0 && i < Nx && j == Ny ) {
		f1xre = (vars[(Ny+1)*(i+1)+j]-vars[(Ny+1)*(i-1)+j])/2*rax;
		f1xim = (vars[(Ny+1)*(i+1)+j+(Nx+1)*(Ny+1)]-vars[(Ny+1)*(i-1)+j+(Nx+1)*(Ny+1)])/2*rax;
		f2xre = (vars[(Ny+1)*(i+1)+j+2*(Nx+1)*(Ny+1)]-vars[(Ny+1)*(i-1)+j+2*(Nx+1)*(Ny+1)])/2*rax;
		f2xim = (vars[(Ny+1)*(i+1)+j+3*(Nx+1)*(Ny+1)]-vars[(Ny+1)*(i-1)+j+3*(Nx+1)*(Ny+1)])/2*rax;

		f1yre = (vars[(Ny+1)*i+j]-vars[(Ny+1)*i+j-1])*ray;
		f1yim = (vars[(Ny+1)*i+j+(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j-1+(Nx+1)*(Ny+1)])*ray;
		f2yre = (vars[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j-1+2*(Nx+1)*(Ny+1)])*ray;
		f2yim = (vars[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j-1+3*(Nx+1)*(Ny+1)])*ray;
	} else if ( i == Nx && j == 0 ) {
		f1xre = (vars[(Ny+1)*i+j]-vars[(Ny+1)*(i-1)+j])*rax;
		f1xim = (vars[(Ny+1)*i+j+(Nx+1)*(Ny+1)]-vars[(Ny+1)*(i-1)+j+(Nx+1)*(Ny+1)])*rax;
		f2xre = (vars[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)]-vars[(Ny+1)*(i-1)+j+2*(Nx+1)*(Ny+1)])*rax;
		f2xim = (vars[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)]-vars[(Ny+1)*(i-1)+j+3*(Nx+1)*(Ny+1)])*rax;

		f1yre = (vars[(Ny+1)*i+j+1]-vars[(Ny+1)*i+j])*ray;
		f1yim = (vars[(Ny+1)*i+j+1+(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j+(Nx+1)*(Ny+1)])*ray;
		f2yre = (vars[(Ny+1)*i+j+1+2*(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)])*ray;
		f2yim = (vars[(Ny+1)*i+j+1+3*(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)])*ray;
	} else if ( i == Nx && j > 0 && j < Ny ) {
		f1xre = (vars[(Ny+1)*i+j]-vars[(Ny+1)*(i-1)+j])*rax;
		f1xim = (vars[(Ny+1)*i+j+(Nx+1)*(Ny+1)]-vars[(Ny+1)*(i-1)+j+(Nx+1)*(Ny+1)])*rax;
		f2xre = (vars[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)]-vars[(Ny+1)*(i-1)+j+2*(Nx+1)*(Ny+1)])*rax;
		f2xim = (vars[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)]-vars[(Ny+1)*(i-1)+j+3*(Nx+1)*(Ny+1)])*rax;

		f1yre = (vars[(Ny+1)*i+j+1]-vars[(Ny+1)*i+j-1])/2*ray;
		f1yim = (vars[(Ny+1)*i+j+1+(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j-1+(Nx+1)*(Ny+1)])/2*ray;
		f2yre = (vars[(Ny+1)*i+j+1+2*(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j-1+2*(Nx+1)*(Ny+1)])/2*ray;
		f2yim = (vars[(Ny+1)*i+j+1+3*(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j-1+3*(Nx+1)*(Ny+1)])/2*ray;
	} else if ( i == Nx && j == Ny ) {
		f1xre = (vars[(Ny+1)*i+j]-vars[(Ny+1)*(i-1)+j])*rax;
		f1xim = (vars[(Ny+1)*i+j+(Nx+1)*(Ny+1)]-vars[(Ny+1)*(i-1)+j+(Nx+1)*(Ny+1)])*rax;
		f2xre = (vars[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)]-vars[(Ny+1)*(i-1)+j+2*(Nx+1)*(Ny+1)])*rax;
		f2xim = (vars[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)]-vars[(Ny+1)*(i-1)+j+3*(Nx+1)*(Ny+1)])*rax;

		f1yre = (vars[(Ny+1)*i+j]-vars[(Ny+1)*i+j-1])*ray;
		f1yim = (vars[(Ny+1)*i+j+(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j-1+(Nx+1)*(Ny+1)])*ray;
		f2yre = (vars[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j-1+2*(Nx+1)*(Ny+1)])*ray;
		f2yim = (vars[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)]-vars[(Ny+1)*i+j-1+3*(Nx+1)*(Ny+1)])*ray;
	}

	double Axf1re = vars[(Ny+1)*i+j+4*(Nx+1)*(Ny+1)]*vars[(Ny+1)*i+j];
	double Axf1im = vars[(Ny+1)*i+j+4*(Nx+1)*(Ny+1)]*vars[(Ny+1)*i+j+(Nx+1)*(Ny+1)];
	double Axf2re = vars[(Ny+1)*i+j+4*(Nx+1)*(Ny+1)]*vars[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)];
	double Axf2im = vars[(Ny+1)*i+j+4*(Nx+1)*(Ny+1)]*vars[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)];
	double Ayf1re = vars[(Ny+1)*i+j+5*(Nx+1)*(Ny+1)]*vars[(Ny+1)*i+j];
	double Ayf1im = vars[(Ny+1)*i+j+5*(Nx+1)*(Ny+1)]*vars[(Ny+1)*i+j+(Nx+1)*(Ny+1)];
	double Ayf2re = vars[(Ny+1)*i+j+5*(Nx+1)*(Ny+1)]*vars[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)];
	double Ayf2im = vars[(Ny+1)*i+j+5*(Nx+1)*(Ny+1)]*vars[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)];

	double T1x = (f1xre-e*Axf1im)*(f1xre-e*Axf1im)+(f1xim+e*Axf1re)*(f1xim+e*Axf1re);
	double T1y = (f1yre-e*Ayf1im)*(f1yre-e*Ayf1im)+(f1yim+e*Ayf1re)*(f1yim+e*Ayf1re);
	double T2x = (f2xre-e*Axf2im)*(f2xre-e*Axf2im)+(f2xim+e*Axf2re)*(f2xim+e*Axf2re);
	double T2y = (f2yre-e*Ayf2im)*(f2yre-e*Ayf2im)+(f2yim+e*Ayf2re)*(f2yim+e*Ayf2re);
	double T = 0.5*((T1x + T1y) + (T2x + T2y));
	return T;
}

// Функция вычисления энергии магнитного поля, связанной с узлом с координатами i, j
inline double CalculateH2( int i, int j, double * vars )
{
	double Ay_up = (vars[(Ny+1)*i+j+5*(Nx+1)*(Ny+1)]+vars[(Ny+1)*i+j+1+5*(Nx+1)*(Ny+1)])/2;
	double Ay_down = (vars[(Ny+1)*(i+1)+j+5*(Nx+1)*(Ny+1)]+vars[(Ny+1)*(i+1)+j+1+5*(Nx+1)*(Ny+1)])/2;
	double Ax_right = (vars[(Ny+1)*i+j+1+4*(Nx+1)*(Ny+1)]+vars[(Ny+1)*(i+1)+j+1+4*(Nx+1)*(Ny+1)])/2;
	double Ax_left = (vars[(Ny+1)*i+j+4*(Nx+1)*(Ny+1)]+vars[(Ny+1)*(i+1)+j+4*(Nx+1)*(Ny+1)])/2;
	double H = (-Ay_up*ay+Ax_left*ax+Ay_down*ay-Ax_right*ax)*rax*ray;

	return 0.5*H*H;
}

// Функция вычисления энергии связанной с узлом с координатами i, j и его соседями
inline double CalculateNeihbourF( int i, int j, double * x )
{
	double F = 0.0;
	for ( int ii = i-1; ii <= i+1; ii++ ) {
		for ( int jj = j-1; jj <= j+1; jj++ ) {
			if ( (ii >= 0) && (ii<Nx+1) && (jj >= 0) && (jj<Ny+1) ) {
				F += CalculateW(ii, jj, x) + CalculateT(ii, jj, x);
			}
		}
	}
	if ( i < Nx && j < Ny && i > 0 && j > 0 ) {
		F += CalculateH2(i-1,j-1,x)+CalculateH2(i,j-1,x)+CalculateH2(i-1,j,x)+CalculateH2(i,j,x);
	} else if ( i == 0 && j == 0 ) {
		F += CalculateH2(i,j, x);
	} else if ( i == 0 && j > 0 && j < Ny ) {
		F += CalculateH2(i,j-1,x)+CalculateH2(i,j,x);
	} else if ( i == 0 && j == Ny ) {
		F += CalculateH2(i,j-1,x);
	} else if ( i > 0 && i < Nx && j == 0 )	{
		F += CalculateH2(i-1,j,x)+CalculateH2(i,j,x);
	} else if ( i > 0 && i < Nx && j == Ny ) {
		F += CalculateH2(i-1,j-1,x)+CalculateH2(i,j-1,x);
	} else if ( i == Nx && j == 0 )	{
		F += CalculateH2(i-1,j,x);
	} else if (i == Nx && j > 0 && j < Ny ) {
		F += CalculateH2(i-1,j-1,x)+CalculateH2(i-1,j,x);
	} else if ( i == Nx && j == Ny ) {
		F += CalculateH2(i-1,j-1,x);
	}
	return ax*ay*F;
}

// Функция вычисления градиента
void CalculateGradient( double * vars, double * grad )
{
	for ( int i = 0; i < Nx+1; i++ ) {
		for ( int j = 0; j < Ny+1; j++ ) {
			//if(!((Nx/2 - i == d)&&(Ny/2 == j) || (Nx/2 - i == - d)&&(Ny/2 == j)))
			{
				//F0 = CalculateNeihbourF(i, j, vars);
				//f1real
				vars[(Ny+1)*i+j] += deltap;
				F1 = CalculateNeihbourF(i, j, vars);
				vars[(Ny+1)*i+j] -= 2*deltap;
				F2 = CalculateNeihbourF(i, j, vars);
				vars[(Ny+1)*i+j] += deltap;
				grad[(Ny+1)*i+j] = 0.5*(F1-F2)*rdeltap;

				//f1imag
				vars[(Ny+1)*i+j+(Nx+1)*(Ny+1)] += deltap;
				F1 = CalculateNeihbourF(i, j, vars);
				vars[(Ny+1)*i+j+(Nx+1)*(Ny+1)] -= 2*deltap;
				F2 = CalculateNeihbourF(i, j, vars);
				vars[(Ny+1)*i+j+(Nx+1)*(Ny+1)] += deltap;
				grad[(Ny+1)*i+j+(Nx+1)*(Ny+1)] = 0.5*(F1-F2)*rdeltap;

				//f2real
				vars[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)] += deltap;
				F1 = CalculateNeihbourF(i, j, vars);
				vars[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)] -= 2*deltap;
				F2 = CalculateNeihbourF(i, j, vars);
				vars[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)] += deltap;
				grad[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)] = 0.5*(F1-F2)*rdeltap;

				//f2imag
				vars[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)] += deltap;
				F1 = CalculateNeihbourF(i, j, vars);
				vars[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)] -= 2*deltap;
				F2 = CalculateNeihbourF(i, j, vars);
				vars[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)] += deltap;
				grad[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)] = 0.5*(F1-F2)*rdeltap;
			}
			//else
			//{
			//	grad[(Ny+1)*i+j] = 0.0;
			//	grad[(Ny+1)*i+j+(Nx+1)*(Ny+1)] = 0.0;
			//	grad[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)] = 0.0;
			//	grad[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)] = 0.0;
			//};

			//Ax
			vars[(Ny+1)*i+j+4*(Nx+1)*(Ny+1)] += deltap;
			F1 = CalculateNeihbourF(i, j, vars);
			vars[(Ny+1)*i+j+4*(Nx+1)*(Ny+1)] -= 2*deltap;
			F2 = CalculateNeihbourF(i, j, vars);
			vars[(Ny+1)*i+j+4*(Nx+1)*(Ny+1)] += deltap;
			grad[(Ny+1)*i+j+4*(Nx+1)*(Ny+1)] = 0.5*(F1-F2)*rdeltap;

			//Ay
			vars[(Ny+1)*i+j+5*(Nx+1)*(Ny+1)] += deltap;
			F1 = CalculateNeihbourF(i, j, vars);
			vars[(Ny+1)*i+j+5*(Nx+1)*(Ny+1)] -= 2*deltap;
			F2 = CalculateNeihbourF(i, j, vars);
			vars[(Ny+1)*i+j+5*(Nx+1)*(Ny+1)] += deltap;
			grad[(Ny+1)*i+j+5*(Nx+1)*(Ny+1)] = 0.5*(F1-F2)*rdeltap;
		}
	}
}

// Заполнение таблиц энергии узлов  и интегрирование методом прямоугольников
inline double FillEnergyTableSquare( double * x )
{
	double F = 0.0;
	double W0 = 0.0;
	double T0 = 0.0;
	double H20 = 0.0;
	for ( int i = 0; i < Nx+1; i++ ) {
		for ( int j = 0; j < Ny+1; j++ ) {
			W[i][j] = CalculateW(i, j, x);
			T[i][j] = CalculateT(i, j, x);
			//W0 += W[i][j];
			//T0 += T[i][j];
			F += T[i][j] + W[i][j];
			//F += CalculateW(i, j) +  CalculateT(i, j);
		}
	}
	for ( int i = 0; i < Nx; i++ ) {
		for ( int j = 0; j < Ny; j++ ) {
			H2[i][j] = CalculateH2(i, j, x);
			//H20 += H2[i][j];
			F += H2[i][j];
			//F += CalculateH2(i, j);
		}
	}
	return ax*ay*F;
}

// Заполнение таблиц энергии узлов  и интегрирование методом трапеций
inline double FillEnergyTableTrapz( double * x )
{
	double F1 = 0.0;
	double F2 = 0.0;
	double W0 = 0.0;
	double T0 = 0.0;
	double H20 = 0.0;
	for ( int i = 0; i < Nx+1; i++ ) {
		for ( int j = 0; j < Ny+1; j++ ) {
			W[i][j] = CalculateW(i, j, x);
			T[i][j] = CalculateT(i, j, x);
			//W0 += W[i][j];
			//T0 += T[i][j];
			//F1 += T[i][j] + W[i][j];
			//F += CalculateW(i, j) +  CalculateT(i, j);
		}
	}
	F1 += T[0][0] + W[0][0];
	F1 += T[Nx][Ny] + W[Nx][Ny];
	F1 += T[0][Ny] + W[0][Ny];
	F1 += T[Nx][0] + W[Nx][0];
	for ( int i = 1; i < Nx; i++ ) {
		F1 += 2*(T[i][0] + W[i][0]);
		F1 += 2*(T[i][Ny] + W[i][Ny]);
	}
	for ( int j = 1; j < Ny; j++ ) {
		F1 += 2*(T[0][j] + W[0][j]);
		F1 += 2*(T[Nx][j] + W[Nx][j]);
	}
	for ( int i = 1; i < Nx; i++ ) {
		for ( int j = 1; j < Ny; j++ ) {
			F1 += 4*(T[i][j] + W[i][j]);
		}
	}
	for ( int i = 0; i < Nx; i++ ) {
		for ( int j = 0; j < Ny; j++ ) {
			H2[i][j] = CalculateH2(i, j, x);
			//H20 += H2[i][j];
			//F2 += H2[i][j];
			//F += CalculateH2(i, j);
		}
	}
	F2 += H2[0][0];
	F2 += H2[Nx-1][Ny-1];
	F2 += H2[0][Ny-1];
	F2 += H2[Nx-1][0];
	for ( int i = 1; i < Nx-1; i++ ) {
		F2 += 2*H2[i][0];
		F2 += 2*H2[i][Ny-1];
	}
	for ( int j = 1; j < Ny-1; j++ ) {
		F2 += 2*H2[0][j];
		F2 += 2*H2[Nx-1][j];
	}
	for ( int i = 1; i < Nx-1; i++ ) {
		for ( int j = 1; j < Ny-1; j++ ) {
			F2 += 4*H2[i][j];
		}
	}
	return 0.25*ax*ay*(F1+F2);
}

// Функция вычисления "потенциальной" энергии, связанной с узлом с координатами i, j
template <class Type>
inline Type CalculateWPT( int i, int j, vector<Type> &varsp )
{
	Type f1real = varsp[(Ny+1)*i+j];
	Type f1imag = varsp[(Ny+1)*i+j+(Nx+1)*(Ny+1)];
	Type f2real = varsp[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)];
	Type f2imag = varsp[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)];
	Type Theta1 = CppAD::atan2(f1imag+1e-30, f1real+1e-30);
	Type Theta2 = CppAD::atan2(f2imag+1e-30, f2real+1e-30);
	Type f1abs = f1real*f1real+f1imag*f1imag + 1e-30;
	Type f2abs = f2real*f2real+f2imag*f2imag + 1e-30;
	Type W = 0;
	W = 0.5*((2*alpha1 + beta1*f1abs)*f1abs + (2*alpha2 + beta2*f2abs)*f2abs) - etta*sqrt(f1abs*f2abs)*(CppAD::cos(Theta2 - Theta1));
	return W;
}

// Функция аналитического вычисления "кинетической" энергии, связанной с узлом с координатами i, j
template <class Type>
inline Type CalculateTPT( int i, int j, vector<Type> &varsp )
{
	Type f1xre;
	Type f2xre;
	Type f1yre;
	Type f2yre;
	Type f1xim;
	Type f2xim;
	Type f1yim;
	Type f2yim;

	if ( i < Nx && j < Ny && i > 0 && j > 0 ) {
		f1xre = (varsp[(Ny+1)*(i+1)+j]-varsp[(Ny+1)*(i-1)+j])/2*rax;
		f1xim = (varsp[(Ny+1)*(i+1)+j+(Nx+1)*(Ny+1)]-varsp[(Ny+1)*(i-1)+j+(Nx+1)*(Ny+1)])/2*rax;
		f2xre = (varsp[(Ny+1)*(i+1)+j+2*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*(i-1)+j+2*(Nx+1)*(Ny+1)])/2*rax;
		f2xim = (varsp[(Ny+1)*(i+1)+j+3*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*(i-1)+j+3*(Nx+1)*(Ny+1)])/2*rax;

		f1yre = (varsp[(Ny+1)*i+j+1]-varsp[(Ny+1)*i+j-1])/2*ray;
		f1yim = (varsp[(Ny+1)*i+j+1+(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j-1+(Nx+1)*(Ny+1)])/2*ray;
		f2yre = (varsp[(Ny+1)*i+j+1+2*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j-1+2*(Nx+1)*(Ny+1)])/2*ray;
		f2yim = (varsp[(Ny+1)*i+j+1+3*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j-1+3*(Nx+1)*(Ny+1)])/2*ray;
	} else if ( i == 0 && j == 0 ) {
		f1xre = (varsp[(Ny+1)*(i+1)+j]-varsp[(Ny+1)*i+j])*rax;
		f1xim = (varsp[(Ny+1)*(i+1)+j+(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j+(Nx+1)*(Ny+1)])*rax;
		f2xre = (varsp[(Ny+1)*(i+1)+j+2*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)])*rax;
		f2xim = (varsp[(Ny+1)*(i+1)+j+3*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)])*rax;

		f1yre = (varsp[(Ny+1)*i+j+1]-varsp[(Ny+1)*i+j])*ray;
		f1yim = (varsp[(Ny+1)*i+j+1+(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j+(Nx+1)*(Ny+1)])*ray;
		f2yre = (varsp[(Ny+1)*i+j+1+2*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)])*ray;
		f2yim = (varsp[(Ny+1)*i+j+1+3*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)])*ray;
	} else if ( i == 0 && j > 0 && j < Ny ) {
		f1xre = (varsp[(Ny+1)*(i+1)+j]-varsp[(Ny+1)*i+j])*rax;
		f1xim = (varsp[(Ny+1)*(i+1)+j+(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j+(Nx+1)*(Ny+1)])*rax;
		f2xre = (varsp[(Ny+1)*(i+1)+j+2*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)])*rax;
		f2xim = (varsp[(Ny+1)*(i+1)+j+3*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)])*rax;

		f1yre = (varsp[(Ny+1)*i+j+1]-varsp[(Ny+1)*i+j-1])/2*ray;
		f1yim = (varsp[(Ny+1)*i+j+1+(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j-1+(Nx+1)*(Ny+1)])/2*ray;
		f2yre = (varsp[(Ny+1)*i+j+1+2*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j-1+2*(Nx+1)*(Ny+1)])/2*ray;
		f2yim = (varsp[(Ny+1)*i+j+1+3*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j-1+3*(Nx+1)*(Ny+1)])/2*ray;
	} else if ( i == 0 && j == Ny ) {
		f1xre = (varsp[(Ny+1)*(i+1)+j]-varsp[(Ny+1)*i+j])*rax;
		f1xim = (varsp[(Ny+1)*(i+1)+j+(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j+(Nx+1)*(Ny+1)])*rax;
		f2xre = (varsp[(Ny+1)*(i+1)+j+2*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)])*rax;
		f2xim = (varsp[(Ny+1)*(i+1)+j+3*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)])*rax;

		f1yre = (varsp[(Ny+1)*i+j]-varsp[(Ny+1)*i+j-1])*ray;
		f1yim = (varsp[(Ny+1)*i+j+(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j-1+(Nx+1)*(Ny+1)])*ray;
		f2yre = (varsp[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j-1+2*(Nx+1)*(Ny+1)])*ray;
		f2yim = (varsp[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j-1+3*(Nx+1)*(Ny+1)])*ray;
	} else if ( i > 0 && i < Nx && j == 0 )	{
		f1xre = (varsp[(Ny+1)*(i+1)+j]-varsp[(Ny+1)*(i-1)+j])/2*rax;
		f1xim = (varsp[(Ny+1)*(i+1)+j+(Nx+1)*(Ny+1)]-varsp[(Ny+1)*(i-1)+j+(Nx+1)*(Ny+1)])/2*rax;
		f2xre = (varsp[(Ny+1)*(i+1)+j+2*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*(i-1)+j+2*(Nx+1)*(Ny+1)])/2*rax;
		f2xim = (varsp[(Ny+1)*(i+1)+j+3*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*(i-1)+j+3*(Nx+1)*(Ny+1)])/2*rax;

		f1yre = (varsp[(Ny+1)*i+j+1]-varsp[(Ny+1)*i+j])*ray;
		f1yim = (varsp[(Ny+1)*i+j+1+(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j+(Nx+1)*(Ny+1)])*ray;
		f2yre = (varsp[(Ny+1)*i+j+1+2*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)])*ray;
		f2yim = (varsp[(Ny+1)*i+j+1+3*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)])*ray;
	} else if ( i > 0 && i < Nx && j == Ny ) {
		f1xre = (varsp[(Ny+1)*(i+1)+j]-varsp[(Ny+1)*(i-1)+j])/2*rax;
		f1xim = (varsp[(Ny+1)*(i+1)+j+(Nx+1)*(Ny+1)]-varsp[(Ny+1)*(i-1)+j+(Nx+1)*(Ny+1)])/2*rax;
		f2xre = (varsp[(Ny+1)*(i+1)+j+2*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*(i-1)+j+2*(Nx+1)*(Ny+1)])/2*rax;
		f2xim = (varsp[(Ny+1)*(i+1)+j+3*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*(i-1)+j+3*(Nx+1)*(Ny+1)])/2*rax;

		f1yre = (varsp[(Ny+1)*i+j]-varsp[(Ny+1)*i+j-1])*ray;
		f1yim = (varsp[(Ny+1)*i+j+(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j-1+(Nx+1)*(Ny+1)])*ray;
		f2yre = (varsp[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j-1+2*(Nx+1)*(Ny+1)])*ray;
		f2yim = (varsp[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j-1+3*(Nx+1)*(Ny+1)])*ray;
	} else if ( i == Nx && j == 0 )	{
		f1xre = (varsp[(Ny+1)*i+j]-varsp[(Ny+1)*(i-1)+j])*rax;
		f1xim = (varsp[(Ny+1)*i+j+(Nx+1)*(Ny+1)]-varsp[(Ny+1)*(i-1)+j+(Nx+1)*(Ny+1)])*rax;
		f2xre = (varsp[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*(i-1)+j+2*(Nx+1)*(Ny+1)])*rax;
		f2xim = (varsp[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*(i-1)+j+3*(Nx+1)*(Ny+1)])*rax;

		f1yre = (varsp[(Ny+1)*i+j+1]-varsp[(Ny+1)*i+j])*ray;
		f1yim = (varsp[(Ny+1)*i+j+1+(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j+(Nx+1)*(Ny+1)])*ray;
		f2yre = (varsp[(Ny+1)*i+j+1+2*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)])*ray;
		f2yim = (varsp[(Ny+1)*i+j+1+3*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)])*ray;
	} else if ( i == Nx && j > 0 && j < Ny ) {
		f1xre = (varsp[(Ny+1)*i+j]-varsp[(Ny+1)*(i-1)+j])*rax;
		f1xim = (varsp[(Ny+1)*i+j+(Nx+1)*(Ny+1)]-varsp[(Ny+1)*(i-1)+j+(Nx+1)*(Ny+1)])*rax;
		f2xre = (varsp[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*(i-1)+j+2*(Nx+1)*(Ny+1)])*rax;
		f2xim = (varsp[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*(i-1)+j+3*(Nx+1)*(Ny+1)])*rax;

		f1yre = (varsp[(Ny+1)*i+j+1]-varsp[(Ny+1)*i+j-1])/2*ray;
		f1yim = (varsp[(Ny+1)*i+j+1+(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j-1+(Nx+1)*(Ny+1)])/2*ray;
		f2yre = (varsp[(Ny+1)*i+j+1+2*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j-1+2*(Nx+1)*(Ny+1)])/2*ray;
		f2yim = (varsp[(Ny+1)*i+j+1+3*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j-1+3*(Nx+1)*(Ny+1)])/2*ray;
	} else if ( i == Nx && j == Ny ) {
		f1xre = (varsp[(Ny+1)*i+j]-varsp[(Ny+1)*(i-1)+j])*rax;
		f1xim = (varsp[(Ny+1)*i+j+(Nx+1)*(Ny+1)]-varsp[(Ny+1)*(i-1)+j+(Nx+1)*(Ny+1)])*rax;
		f2xre = (varsp[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*(i-1)+j+2*(Nx+1)*(Ny+1)])*rax;
		f2xim = (varsp[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*(i-1)+j+3*(Nx+1)*(Ny+1)])*rax;

		f1yre = (varsp[(Ny+1)*i+j]-varsp[(Ny+1)*i+j-1])*ray;
		f1yim = (varsp[(Ny+1)*i+j+(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j-1+(Nx+1)*(Ny+1)])*ray;
		f2yre = (varsp[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j-1+2*(Nx+1)*(Ny+1)])*ray;
		f2yim = (varsp[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)]-varsp[(Ny+1)*i+j-1+3*(Nx+1)*(Ny+1)])*ray;
	}

	Type Axf1re = varsp[(Ny+1)*i+j+4*(Nx+1)*(Ny+1)]*varsp[(Ny+1)*i+j];
	Type Axf1im = varsp[(Ny+1)*i+j+4*(Nx+1)*(Ny+1)]*varsp[(Ny+1)*i+j+(Nx+1)*(Ny+1)];
	Type Axf2re = varsp[(Ny+1)*i+j+4*(Nx+1)*(Ny+1)]*varsp[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)];
	Type Axf2im = varsp[(Ny+1)*i+j+4*(Nx+1)*(Ny+1)]*varsp[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)];
	Type Ayf1re = varsp[(Ny+1)*i+j+5*(Nx+1)*(Ny+1)]*varsp[(Ny+1)*i+j];
	Type Ayf1im = varsp[(Ny+1)*i+j+5*(Nx+1)*(Ny+1)]*varsp[(Ny+1)*i+j+(Nx+1)*(Ny+1)];
	Type Ayf2re = varsp[(Ny+1)*i+j+5*(Nx+1)*(Ny+1)]*varsp[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)];
	Type Ayf2im = varsp[(Ny+1)*i+j+5*(Nx+1)*(Ny+1)]*varsp[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)];

	Type T1x = (f1xre-e*Axf1im)*(f1xre-e*Axf1im)+(f1xim+e*Axf1re)*(f1xim+e*Axf1re);
	Type T1y = (f1yre-e*Ayf1im)*(f1yre-e*Ayf1im)+(f1yim+e*Ayf1re)*(f1yim+e*Ayf1re);
	Type T2x = (f2xre-e*Axf2im)*(f2xre-e*Axf2im)+(f2xim+e*Axf2re)*(f2xim+e*Axf2re);
	Type T2y = (f2yre-e*Ayf2im)*(f2yre-e*Ayf2im)+(f2yim+e*Ayf2re)*(f2yim+e*Ayf2re);
	Type T = 0.5*((T1x + T1y) + (T2x + T2y));
	return T;
}

// Функция аналитического вычисления энергии магнитного поля, связанной с узлом с координатами i, j
template <class Type>
inline Type CalculateH2PT( int i, int j, vector<Type> &varsp )
{
	Type H;
	Type Ay_up = (varsp[(Ny+1)*i+j+5*(Nx+1)*(Ny+1)]+varsp[(Ny+1)*i+j+1+5*(Nx+1)*(Ny+1)])/2;
	Type Ay_down = (varsp[(Ny+1)*(i+1)+j+5*(Nx+1)*(Ny+1)]+varsp[(Ny+1)*(i+1)+j+1+5*(Nx+1)*(Ny+1)])/2;
	Type Ax_right = (varsp[(Ny+1)*i+j+1+4*(Nx+1)*(Ny+1)]+varsp[(Ny+1)*(i+1)+j+1+4*(Nx+1)*(Ny+1)])/2;
	Type Ax_left = (varsp[(Ny+1)*i+j+4*(Nx+1)*(Ny+1)]+varsp[(Ny+1)*(i+1)+j+4*(Nx+1)*(Ny+1)])/2;
	H = (-Ay_up*ay+Ax_left*ax+Ay_down*ay-Ax_right*ax)*rax*ray;

	return 0.5*H*H;
}

// Аналитическое интегрирование методом прямоугольников
template <class Type>
inline Type EnergySquareT( vector<Type> &varsp )
{
	Type F = 0.0;
	Type W0 = 0.0;
	Type T0 = 0.0;
	Type H20 = 0.0;

	for ( int i = 0; i < Nx + 1; i++ ) {
		for ( int j = 0; j < Ny + 1; j++ ) {
			F += CalculateWPT(i, j, varsp) + CalculateTPT(i, j, varsp);
		}
	}
	for ( int i = 0; i < Nx; i++ ) {
		for ( int j = 0; j < Ny; j++ ) {
			F += CalculateH2PT(i, j, varsp);
		}
	}
	return ax*ay*F;
}

// Аналитическое интегрирование методом трапеций
template <class Type>
inline Type EnergyTrapzT(vector<Type> &varsp)
{
	Type F1 = 0.0;
	Type F2 = 0.0;
	Type W0 = 0.0;
	Type T0 = 0.0;
	Type H20 = 0.0;

	F1 += CalculateWP(0, 0, varsp) + CalculateTP(0, 0, varsp);
	F1 += CalculateWP(Nx, Ny, varsp) + CalculateTP(Nx, Ny, varsp);
	F1 += CalculateWP(0, Ny, varsp) + CalculateTP(0, Ny, varsp);
	F1 += CalculateWP(Nx, 0, varsp) + CalculateTP(Nx, 0, varsp);
	for ( int i = 1; i < Nx; i++ ) {
		F1 += 2*(CalculateWP(i, 0, varsp) + CalculateTP(i, 0, varsp));
		F1 += 2*(CalculateWP(i, Ny, varsp) + CalculateTP(i, Ny, varsp));
	}
	for ( int j = 1; j < Ny; j++ ) {
		F1 += 2*(CalculateWP(0, j, varsp) + CalculateTP(0, j, varsp));
		F1 += 2*(CalculateWP(Nx, j, varsp) + CalculateTP(Nx, j, varsp));
	}
	for ( int i = 1; i < Nx; i++ ) {
		for ( int j = 1; j < Ny; j++ ) {
			F1 += 4*(CalculateWP(i, j, varsp) + CalculateTP(i, j, varsp));
		}
	}
	F2 += CalculateH2P(0, 0, varsp);
	F2 += CalculateH2P(Nx-1, Ny-1, varsp);
	F2 += CalculateH2P(0, Ny-1, varsp);
	F2 += CalculateH2P(Nx-1, 0, varsp);
	for ( int i = 1; i < Nx - 1; i++ ) {
		F2 += 2*CalculateH2P(i, 0, varsp);
		F2 += 2*CalculateH2P(i, Ny-1, varsp);
	}
	for ( int j = 1; j < Ny - 1; j++ ) {
		F2 += 2*CalculateH2P(0, j, varsp);
		F2 += 2*CalculateH2P(Nx-1, j, varsp);
	}
	for ( int i = 1; i < Nx - 1; i++ ) {
		for ( int j = 1; j < Ny - 1; j++ ) {
			F2 += 4*CalculateH2P(i, j, varsp);
		}
	}
	return 0.25*ax*ay*(F1+F2);
}

// Блок для минимизации алгоритмом LBFGS из HLBFGS
void evalfunc( int n, double* x, double *prev_x, double* f, double* g )
{
	double ff = FillEnergyTableSquare(x);
	*f = ff;
	vector<double> g1(n);
	vector<double> xx(n);
	for ( int i = 0; i < n; i++ ) {
		xx[i] = x[i];
	}
	g1 = fun.Jacobian( xx );
	for ( int i = 0; i < n; i++ ) {
		g[i] = g1[i];
	}
	g1.clear();
}

void newiteration(int iter, int call_iter, double *x, double* f, double *g,  double* gnorm)
{
	nextTime = get_ticks();
	printf( "[%u][%d]: %.8f %.8f %.8f\n", nextTime - predTime, call_iter, *g, *f, *gnorm );
	predTime = nextTime;
}

int main(int argc, char** argv)
{
	// Расстояние от вихря до начала координат
	p_block_t block[] = {
		{ "alpha1", T_DOUBLE, 0, 0, &alpha1 },
		{ "alpha2", T_DOUBLE, 0, 0, &alpha2 },
		{ "ax", T_DOUBLE, 0, 0, &ax },
		{ "ay", T_DOUBLE, 0, 0, &ay },
		{ "beta1", T_DOUBLE, 0, 0, &beta1 },
		{ "beta2", T_DOUBLE, 0, 0, &beta2 },
		{ "d", T_INT, &d },
		{ "deltap", T_DOUBLE, 0, 0, &deltap },
		{ "e", T_DOUBLE, 0, 0, &e },
		{ "epsilon", T_DOUBLE, 0, 0, &epsilon },
		{ "etta", T_DOUBLE, 0, 0, &etta },
		{ "Nx", T_INT, &Nx },
		{ "Ny", T_INT, &Ny },
		{ NULL }
	};
	int k_file = 0;
	vector<string> names;
	string word;
	ifstream file;

	file.open( argv[1] );
	while ( file >> word ) {
		cout << word << endl;
		names.push_back( word );
	}
    file.close();
    printf( ">> Start working [names = %lu]\n", names.size() );
	while ( names.size() > (size_t) k_file ) {
        printf( ">> Start Job#%02d [ %s ]\n", k_file, names[k_file].c_str() );
        printf( ">> Time start: " );
		get_time();
		config_parser( names[k_file].c_str(), 13, block );
		rax = 1/ax;
		ray = 1/ay;
		iv1 = Nx/2 - d;
		jv1 = Ny/2;
		iv2 = Nx/2 + d;
		jv2 = Ny/2;
		// Создание массивов

		n = 6*(Nx+1)*(Ny+1);
		X.resize( n );
		x.setlength( n );
		df1re = new double* [Nx + 1];
		df2re = new double* [Nx + 1];
		df1im = new double* [Nx + 1];
		df2im = new double* [Nx + 1];
		W = new double* [Nx + 1];
		T = new double* [Nx + 1];
		vars = new double [n];
		grad = new double [n];
		H2 = new double* [Nx +1];
		dAx = new double* [Nx + 1];
		dAy = new double* [Nx + 1];
		for ( int i = 0; i < Nx + 1; i++ ) {
			df1re[i] = new double [Ny + 1];
			df2re[i] = new double [Ny + 1];
			df1im[i] = new double [Ny + 1];
			df2im[i] = new double [Ny + 1];
			W[i] = new double [Ny + 1];
			T[i] = new double [Ny + 1];
			dAx[i] = new double [Ny + 1];
			dAy[i] = new double [Ny + 1];
			H2[i] = new double [Ny + 1];
		}
		for ( int i = 0; i < Nx + 1; i++ ) {
			for ( int j = 0; j < Ny + 1; j++ ) {
				z1 = (Nx / 2 - i - d) * ax + I * ay * (double)(Ny / 2 - j);
				z2 = (Nx / 2 - i + d) * ax + I * ay * (double)(Ny / 2 - j);
				theta1 = arg( z1 );
				theta2 = arg( z2 );
				r1 = abs( z1 );
				r2 = abs( z2 );
				// Начальная конфигурация параметров порядка
				if (z1 == 0.0 || z2 == 0.0) {
					vars[(Ny + 1) * i + j] = 0.0;
					vars[(Ny + 1) * i + j + (Nx + 1) * (Ny + 1)] = 0.0;
					vars[(Ny + 1) * i + j + 2 * (Nx + 1) * (Ny + 1)] = 0.0;
					vars[(Ny + 1) * i + j + 3 * (Nx + 1) * (Ny + 1)] = 0.0;
				} else {
					f1 = u10 * sqrt(0.5 + 0.5 * tanh(4.0 / ksi1 * (r1 - ksi1))) *
						exp(I * theta1) * sqrt(0.5 + 0.5 * tanh(4.0 / ksi1 * (r2 - ksi1))) *
						exp(I * theta2);
					f2 = u20 * sqrt(0.5 + 0.5 * tanh(4.0 / ksi2 * (r2 - ksi2))) *
						exp(I * theta2) * sqrt(0.5 + 0.5 * tanh(4.0 / ksi2 * (r1 - ksi2))) *
						exp(I * theta1);
					vars[(Ny + 1) * i + j] = f1.real();
					vars[(Ny + 1) * i + j + (Nx + 1) * (Ny + 1)] = f1.imag();
					vars[(Ny + 1) * i + j + 2 * (Nx + 1) * (Ny + 1)] = f2.real();
					vars[(Ny + 1) * i + j + 3 * (Nx + 1) * (Ny + 1)] = f2.imag();
					bb = vars[Ny * 0 + 0 + (Nx + 1)*(Ny + 1)];
				}
				// Начальная конфигурация вектор-потенциала Ax
				vars[(Ny + 1) * i + j + 4 * (Nx + 1) * (Ny + 1)] = 1.0 / e / (r1 + r2) *
					sin(theta1 + theta2);
				// Начальная конфигурация вектор-потенциала Ay
				vars[(Ny + 1) * i + j + 5 * (Nx + 1) * (Ny + 1)] = -1.0 / e / (r1 + r2) *
					cos(theta1 + theta2);
			}
		}
		//Конструирование функции для автоматического дифференцирования
		for ( int i = 0; i < n; i++ ) {
			X[i] = 1.0;
		}
		CppAD::Independent( X );
		size_t m = 1;
		Y.resize( m );
		Y[0] = EnergySquareT( X );
		fun = CppAD::ADFun<double>( X, Y );
		
		//Конец конструирования функции для автоматического дифференцирования
        printf( ">> Start energy is %+.8f\n", FillEnergyTableSquare( vars ) );
		double parameter[20];
		parameter[6] = 0.01;
		int info[20];
		// initialize
		INIT_HLBFGS( parameter, info );
		// info[4] = num_iter;
		info[6] = 0;
		info[7] = 0;
		info[10] = 0;
		info[11] = 1;
        nextTime = get_ticks();
		HLBFGS( n, 5, vars, evalfunc, 0, HLBFGS_UPDATE_Hessian, newiteration, parameter, info );
		// Сохранение результатов
		// en1 = FillEnergyTableSquare(vars);

		char buffer[32];
		sprintf( buffer, "job#%02d_%s", k_file, names[k_file].c_str() );
#ifdef WIN32
		mkdir( buffer );
#else
        mkdir( buffer, S_IRWXU | S_IRGRP | S_IROTH );
#endif
		chdir( buffer );

        FILE *f;
        f = fopen( "dataF1r.dat", "w" );
        if ( f == NULL ) {
            f = stdout;
            fprintf( f, "\n--- dataF1r ---\n\n" );
        }
		for ( int i = 0; i < Nx + 1; i++ ) {
			for ( int j = 0; j < Ny + 1; j++ ) {
                fprintf( f, "%+.16f ", vars[(Ny+1)*i+j] );
			}
            fprintf( f, "\n" );
		}
        if ( f != stdout ) {
            fclose( f );
        }
		f = fopen( "dataF1c.dat", "w" );
        if ( f == NULL ) {
            f = stdout;
            fprintf( f, "\n--- dataF1c ---\n\n" );
        }
		for ( int i = 0; i < Nx + 1; i++ ) {
			for ( int j = 0; j < Ny + 1; j++ ) {
				fprintf( f, "%+.16f ", vars[(Ny+1)*i+j+(Nx+1)*(Ny+1)] );
			}
			fprintf( f, "\n" );
		}
		if ( f != stdout ) {
            fclose( f );
        }
		f = fopen( "dataF2r.dat", "w" );
        if ( f == NULL ) {
            f = stdout;
            fprintf( f, "\n--- dataF2r ---\n\n" );
        }
		for ( int i = 0; i < Nx + 1; i++ ) {
			for ( int j = 0; j < Ny + 1; j++) {
				fprintf( f, "%+.16f ", vars[(Ny+1)*i+j+2*(Nx+1)*(Ny+1)] );
			}
			fprintf( f, "\n" );
		}
		if ( f != stdout ) {
            fclose( f );
        }
		f = fopen( "dataF2c.dat", "w" );
		if ( f == NULL ) {
            f = stdout;
            fprintf( f, "\n--- dataF2c ---\n\n" );
        }
		for ( int i = 0; i < Nx + 1; i++ ) {
			for ( int j = 0; j < Ny + 1; j++) {
				fprintf( f, "%+.16f ",vars[(Ny+1)*i+j+3*(Nx+1)*(Ny+1)] );
			}
			fprintf( f, "\n" );
		}
		if ( f != stdout ) {
            fclose( f );
        }
		f = fopen( "dataAx.dat", "w" );
		if ( f == NULL ) {
            f = stdout;
            fprintf( f, "\n--- dataAx ---\n\n" );
        }
		for ( int i = 0; i < Nx + 1; i++ ) {
			for ( int j = 0; j < Ny + 1; j++) {
				fprintf( f, "%+.16f ", vars[(Ny+1)*i+j+4*(Nx+1)*(Ny+1)] );
			}
			fprintf( f, "\n" );
		}
		if ( f != stdout ) {
            fclose( f );
        }
		f = fopen( "dataAy.dat", "w" );
		if ( f == NULL ) {
            f = stdout;
            fprintf( f, "\n--- dataAy ---\n\n" );
        }		
        for ( int i = 0; i < Nx + 1; i++ ) {
			for ( int j = 0; j < Ny + 1; j++ ) {
				fprintf( f, "%+.16f ", vars[(Ny+1)*i+j+5*(Nx+1)*(Ny+1)] );
			}
			fprintf( f, "\n" );
		}
		if ( f != stdout ) {
            fclose( f );
        }
		f = fopen( "dataB.dat", "w" );
		if ( f == NULL ) {
            f = stdout;
            fprintf( f, "\n--- dataB ---\n\n" );
        }
		for ( int i = 0; i < Nx + 1; i++ ) {
			for ( int j = 0; j < Ny + 1; j++ ) {
				fprintf( f, "%+.16f ", sqrt( 2 * H2[i][j] ) );
			}
			fprintf( f, "\n" );
		}
		if ( f != stdout ) {
            fclose( f );
        }
		f = fopen( "dataW.dat", "w" );
		if ( f == NULL ) {
            f = stdout;
            fprintf( f, "\n--- dataW ---\n\n" );
        }
		for ( int i = 0; i < Nx + 1; i++ ) {
			for ( int j = 0; j < Ny + 1; j++ ) {
				fprintf( f, "%+.16f ", W[i][j] );
			}
			fprintf( f, "\n" );
		}
		if ( f != stdout ) {
            fclose( f );
        }
		f = fopen( "dataT.dat", "w" );
		if ( f == NULL ) {
            f = stdout;
            fprintf( f, "\n--- dataT ---\n\n" );
        }
		for ( int i = 0; i < Nx + 1; i++ ) {
			for ( int j = 0; j < Ny + 1; j++ ) {
				fprintf( f, "%+.16f ", T[i][j] );
			}
			fprintf( f, "\n" );
		}
		if ( f != stdout ) {
            fclose( f );
        }
		chdir( ".." );
		k_file++;
		k = n = 0;
		for ( int i = 0; i < Nx + 1; i++ ) {
			free( df1re[i] );
			free( df2re[i] );
			free( df1im[i] );
			free( df2im[i] );
			free( W[i] );
			free( T[i] );
			free( dAx[i] );
			free( dAy[i] );
			free( H2[i] );
		}
		free( vars );
		free( grad );
		free( df1re );
		free( df2re );
		free( df1im );
		free( df2im );
		free( W );
		free( T );
		free( H2 );
		free( dAx );
		free( dAy );
        printf( ">> Time end: " );
		get_time();
	}
	return 0;
}
