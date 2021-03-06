/* $Id: testvector.hpp 2576 2012-11-17 13:44:48Z bradbell $ */
# ifndef CPPAD_TESTVECTOR_INCLUDED
# define CPPAD_TESTVECTOR_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-12 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    GNU General Public License Version 3.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
$begin testvector$$
$spell
	CppAD
	cmake
	testvector
	cppad
	Eigen
	ifdef
	hpp
	std
	endif
	ublas
$$

$index CPPAD_TESTVECTOR$$
$index vector, test$$
$index test, vector$$

$section Using The CppAD Test Vector Template Class$$

$head Syntax$$
$codei%CPPAD_TESTVECTOR(%Scalar%)
%$$

$head Introduction$$
Many of the CppAD $cref/examples/example/$$ and tests use 
the $code CPPAD_TESTVECTOR$$ template class to pass information to CppAD.
This is not a true template class because it's syntax uses
$codei%(%Scalar%)%$$ instead of $codei%<%Scalar%>%$$.
This enables us to use
$codei%
	Eigen::Matrix<%Scalar%, Eigen::Dynamic, 1>
%$$
as one of the possible cases for this 'template class'.

$head CppAD::vector$$
If in the $cref/cmake command/cmake/CMake Command/$$
you specify $cref cppad_testvector$$ to be $code cppad$$,
$code CPPAD_CPPADVECTOR$$ will be true.
In this case,
$code CPPAD_TESTVECTOR$$ is defined by the following source code:
$codep */
# if CPPAD_CPPADVECTOR
# define CPPAD_TESTVECTOR(Scalar) CppAD::vector< Scalar >
# endif
/* $$
In this case CppAD will use its own vector for 
many of its examples and tests.

$head std::vector$$
If in the cmake command
you specify $icode cppad_testvector$$ to be $code std$$,
$code CPPAD_STDVECTOR$$ will be true.
In this case,
$code CPPAD_TESTVECTOR$$ is defined by the following source code:
$codep */
# if CPPAD_STDVECTOR
# include <vector>
# define CPPAD_TESTVECTOR(Scalar) std::vector< Scalar >
# endif
/* $$
In this case CppAD will use standard vector for 
many of its examples and tests.

$head boost::numeric::ublas::vector$$
If in the cmake command
you specify $icode cppad_testvector$$ to be $code boost$$,
$code CPPAD_BOOSTVECTOR$$ will be true.
In this case,
$code CPPAD_TESTVECTOR$$ is defined by the following source code:
$codep */
# if CPPAD_BOOSTVECTOR
# include <boost/numeric/ublas/vector.hpp>
# define CPPAD_TESTVECTOR(Scalar) boost::numeric::ublas::vector< Scalar >
# endif
/* $$
In this case CppAD will use this boost vector for 
many of its examples and tests.

$head Eigen Vectors$$
If in the cmake command
you specify $icode cppad_testvector$$ to be $code eigen$$,
$code CPPAD_EIGENVECTOR$$ will be true.
In this case,
$code CPPAD_TESTVECTOR$$ is defined by the following source code:
$codep */
# if CPPAD_EIGENVECTOR
# include <cppad/example/cppad_eigen.hpp>
# define CPPAD_TESTVECTOR(Scalar) Eigen::Matrix< Scalar , Eigen::Dynamic, 1>
# endif
/* $$
In this case CppAD will use the Eigen vector 
for many of its examples and tests.

$end
------------------------------------------------------------------------ 
*/

# endif
