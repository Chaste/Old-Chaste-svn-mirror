/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


#ifndef UBLASCUSTOMFUNCTIONS_HPP_
#define UBLASCUSTOMFUNCTIONS_HPP_

/**
 * @file
 * A collection of useful functions extending the functionality of the
 * Boost Ublas library.
 */

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>

#include "petsc.h"
#include "petscblaslapack.h"
//Promote universal LAPACK name if it's an old version of PETSc
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
#define LAPACKgeev_ LAgeev_
#endif

#include <cfloat>
#define TINY DBL_EPSILON

#include "Exception.hpp"
#include "PetscTools.hpp"



using namespace boost::numeric::ublas;
/**
 * Get the determinant of a ublas matrix
 * N.B. Do not call with non-square matrix
 *
 *
 * @param rM The matrix of which to find the determinant. Must be square.
 */
template<class T>
T Determinant(const boost::numeric::ublas::c_matrix<T, 1, 1>& rM)
{
    using namespace boost::numeric::ublas;

    return rM(0,0);
}
/**
 * Return the determinant of a submatrix after removing a particular row and column
 */
template<class T>
T SubDeterminant(const boost::numeric::ublas::c_matrix<T, 1, 1>& rM, const unsigned missrow, const unsigned misscol)
{
    using namespace boost::numeric::ublas;

    assert(missrow==0);
    assert(misscol==0);
    return 1.0;
}

template<class T>
T Determinant(const boost::numeric::ublas::c_matrix<T,2,2>& rM)
{
    using namespace boost::numeric::ublas;
#define SMALL_POW(x, a)
#define SMALL_POW(x, a)
#define SMALL_POW(x, a)

    return rM(0,0)*rM(1,1) - rM(1,0)*rM(0,1);
}

/**
 * Return the determinant of a submatrix after removing a particular row and column
 */
template<class T>
T SubDeterminant(const boost::numeric::ublas::c_matrix<T, 2, 2>& rM, const unsigned missrow, const unsigned misscol)
{
    using namespace boost::numeric::ublas;

    assert(missrow < 2);
    assert(misscol < 2);

    unsigned row = (missrow==1) ? 0 : 1;
    unsigned col = (misscol==1) ? 0 : 1;
    return rM(row,col);
}

template<class T>
T SubDeterminant(const boost::numeric::ublas::c_matrix<T, 2, 1>& rM, const unsigned missrow, const unsigned misscol)
{
    using namespace boost::numeric::ublas;

    assert(missrow < 2);
    assert(misscol < 1);

    unsigned row = (missrow==1) ? 0 : 1;
    unsigned col = (misscol==1) ? 0 : 1;
    return rM(row,col);
}

template<class T>
T Determinant(const boost::numeric::ublas::c_matrix<T, 3, 3>& rM)
{
    using namespace boost::numeric::ublas;

    return    rM(0,0) * (rM(1,1)*rM(2,2) - rM(1,2)*rM(2,1))
            - rM(0,1) * (rM(1,0)*rM(2,2) - rM(1,2)*rM(2,0))
            + rM(0,2) * (rM(1,0)*rM(2,1) - rM(1,1)*rM(2,0));
}

/**
 * Calculate the generalized determinant of a 2x1 matrix.
 *
 * The generalized determinant is given by det(T) = sqrt(det(T'T));
 */
template<class T>
T Determinant(const boost::numeric::ublas::c_matrix<T, 2, 1>& rM)
{
    using namespace boost::numeric::ublas;

    return   std::sqrt(rM(0,0) * rM(0,0) + rM(1,0) * rM(1,0));
}

/**
 * Calculate the generalized determinant of a 3x1 matrix.
 *
 * The generalized determinant is given by det(T) = sqrt(det(T'T));
 */
template<class T>
T Determinant(const boost::numeric::ublas::c_matrix<T, 3, 1>& rM)
{
    using namespace boost::numeric::ublas;

    return std::sqrt(rM(0,0)*rM(0,0) + rM(1,0)*rM(1,0) + rM(2,0)*rM(2,0));
}

/**
 * Calculate the generalized determinant of a 3x2 matrix.
 *
 * The generalized determinant is given by det(T) = sqrt(det(T'T));
 */
template<class T>
T Determinant(const boost::numeric::ublas::c_matrix<T, 3, 2>& rM)
{
    using namespace boost::numeric::ublas;

    c_matrix<T,2,2> product = prod(trans(rM), rM);

    return std::sqrt(Determinant(product));
}

/**
 * Calculate the generalized determinant of a 3x0 matrix.
 *
 * The generalized determinant is given by det(T) = sqrt(det(T'T));
 *
 * \todo not implemented yet
 */
template<class T>
T Determinant(const boost::numeric::ublas::c_matrix<T, 3, 0>& rM)
{
    NEVER_REACHED;
}

/**
 * Calculate the generalized determinant of a 2x0 matrix.
 *
 * The generalized determinant is given by det(T) = sqrt(det(T'T));
 *
 * \todo not implemented yet
 */
template<class T>
T Determinant(const boost::numeric::ublas::c_matrix<T, 2, 0>& rM)
{
    NEVER_REACHED;
}

/**
 * Calculate the generalized determinant of a 1x0 matrix.
 *
 * The generalized determinant is given by det(T) = sqrt(det(T'T));
 *
 * \todo not implemented yet
 */
template<class T>
T Determinant(const boost::numeric::ublas::c_matrix<T, 1, 0>& rM)
{
    NEVER_REACHED;
}

#if defined(__xlC__)
/* IBM compiler doesn't support zero-sized arrays*/
#else
/**
 * Return the determinant of a submatrix after removing a particular row and column
 */
template<class T>
T SubDeterminant(const boost::numeric::ublas::c_matrix<T, 3, 0>& rM, const unsigned missrow, const unsigned misscol)
{
    NEVER_REACHED;
}

/**
 * Return the determinant of a submatrix after removing a particular row and column
 */
template<class T>
T SubDeterminant(const boost::numeric::ublas::c_matrix<T, 2, 0>& rM, const unsigned missrow, const unsigned misscol)
{
    NEVER_REACHED;
}

/**
 * Return the determinant of a submatrix after removing a particular row and column
 */
template<class T>
T SubDeterminant(const boost::numeric::ublas::c_matrix<T, 1, 0>& rM, const unsigned missrow, const unsigned misscol)
{
    NEVER_REACHED;
}
#endif

/**
 * Return the determinant of a submatrix after removing a particular row and column
 */
template<class T>
T SubDeterminant(const boost::numeric::ublas::c_matrix<T, 3, 2>& rM, const unsigned missrow, const unsigned misscol)
{
    using namespace boost::numeric::ublas;

    assert(missrow < 3);
    //assert(misscol < 2);

    unsigned lorow = (missrow==0) ? 1 : 0;
    unsigned hirow = (missrow==2) ? 1 : 2;
    unsigned locol = (misscol==0) ? 1 : 0;
    unsigned hicol = (misscol==2) ? 1 : 2;
    return rM(lorow,locol)*rM(hirow,hicol) - rM(lorow,hicol)*rM(hirow,locol);
}

/**
 * Return the determinant of a submatrix after removing a particular row and column
 */
template<class T>
T SubDeterminant(const boost::numeric::ublas::c_matrix<T, 3, 1>& rM, const unsigned missrow, const unsigned misscol)
{
    using namespace boost::numeric::ublas;

    assert(missrow < 3);
    assert(misscol < 1);

    unsigned lorow = (missrow==0) ? 1 : 0;
    unsigned hirow = (missrow==2) ? 1 : 2;
    unsigned locol = (misscol==0) ? 1 : 0;
    unsigned hicol = (misscol==2) ? 1 : 2;
    return rM(lorow,locol)*rM(hirow,hicol) - rM(lorow,hicol)*rM(hirow,locol);
}

/**
 * Return the determinant of a submatrix after removing a particular row and column
 */
template<class T>
T SubDeterminant(const boost::numeric::ublas::c_matrix<T, 3, 3>& rM, const unsigned missrow, const unsigned misscol)
{
    using namespace boost::numeric::ublas;

    assert(missrow < 3);
    assert(misscol < 3);

    unsigned lorow = (missrow==0) ? 1 : 0;
    unsigned hirow = (missrow==2) ? 1 : 2;
    unsigned locol = (misscol==0) ? 1 : 0;
    unsigned hicol = (misscol==2) ? 1 : 2;
    return rM(lorow,locol)*rM(hirow,hicol) - rM(lorow,hicol)*rM(hirow,locol);
}

/**
 * Get the inverse of a ublas matrix
 * N.B. Do not call with non-square matrix
 *
 * @param m The matrix of which to find the inverse. Must be square.
 */
template<class T>
boost::numeric::ublas::c_matrix<T, 1, 1> Inverse(const boost::numeric::ublas::c_matrix<T, 1, 1>& rM)
{
    using namespace boost::numeric::ublas;

    c_matrix<T,1,1> inverse;
    T det = Determinant(rM);
    assert( fabs(det) > TINY ); // else it is a singular matrix
    inverse(0,0) =  1.0/det;
    return inverse;
}

template<class T>
boost::numeric::ublas::c_matrix<T, 2, 3> Inverse(const boost::numeric::ublas::c_matrix<T, 3, 2>& rM)
{
    using namespace boost::numeric::ublas;

    c_matrix<T, 2, 3> inverse;

    //
    // calculate (T'T)^-1, where T'T = (a b)
    //                                 (c d)

    T a = rM(0,0)*rM(0,0) + rM(1,0)*rM(1,0) + rM(2,0)*rM(2,0);
    T b = rM(0,0)*rM(0,1) + rM(1,0)*rM(1,1) + rM(2,0)*rM(2,1);
    T c = b;
    T d = rM(0,1)*rM(0,1) + rM(1,1)*rM(1,1) + rM(2,1)*rM(2,1);

    T det = a*d - b*c;

    T a_inv =  d/det;
    T b_inv = -b/det;
    T c_inv = -c/det;
    T d_inv =  a/det;

    inverse(0,0) = a_inv*rM(0,0) + b_inv*rM(0,1);
    inverse(1,0) = c_inv*rM(0,0) + d_inv*rM(0,1);
    inverse(0,1) = a_inv*rM(1,0) + b_inv*rM(1,1);
    inverse(1,1) = c_inv*rM(1,0) + d_inv*rM(1,1);
    inverse(0,2) = a_inv*rM(2,0) + b_inv*rM(2,1);
    inverse(1,2) = c_inv*rM(2,0) + d_inv*rM(2,1);

    return inverse;
}

template<class T>
boost::numeric::ublas::c_matrix<T, 2, 2> Inverse(const boost::numeric::ublas::c_matrix<T, 2, 2>& rM)
{
    using namespace boost::numeric::ublas;

    c_matrix<T, 2, 2> inverse;
    T det = Determinant(rM);

    assert( fabs(det) > TINY ); // else it is a singular matrix
    inverse(0,0)  =  rM(1,1)/det;
    inverse(0,1)  = -rM(0,1)/det;
    inverse(1,0)  = -rM(1,0)/det;
    inverse(1,1)  =  rM(0,0)/det;
    return inverse;
}

template<class T>
boost::numeric::ublas::c_matrix<T, 3, 3> Inverse(const boost::numeric::ublas::c_matrix<T, 3, 3>& rM)
{
    using namespace boost::numeric::ublas;

    c_matrix<T, 3, 3> inverse;
    T det = Determinant(rM);
    assert( fabs(det) > TINY ); // else it is a singular matrix

    inverse(0,0) =  (rM(1,1)*rM(2,2) - rM(1,2)*rM(2,1))/det;
    inverse(1,0) = -(rM(1,0)*rM(2,2) - rM(1,2)*rM(2,0))/det;
    inverse(2,0) =  (rM(1,0)*rM(2,1) - rM(1,1)*rM(2,0))/det;
    inverse(0,1) = -(rM(0,1)*rM(2,2) - rM(0,2)*rM(2,1))/det;
    inverse(1,1) =  (rM(0,0)*rM(2,2) - rM(0,2)*rM(2,0))/det;
    inverse(2,1) = -(rM(0,0)*rM(2,1) - rM(0,1)*rM(2,0))/det;
    inverse(0,2) =  (rM(0,1)*rM(1,2) - rM(0,2)*rM(1,1))/det;
    inverse(1,2) = -(rM(0,0)*rM(1,2) - rM(0,2)*rM(1,0))/det;
    inverse(2,2) =  (rM(0,0)*rM(1,1) - rM(0,1)*rM(1,0))/det;

    return inverse;
}

/**
 * Calculates the pseudo-inverse of a 2x1 matrix.
 *
 * The pseudo inverse is given by pinv(T) = (T'T)^(-1)*T'
 */
template<class T>
boost::numeric::ublas::c_matrix<T, 1, 2> Inverse(const boost::numeric::ublas::c_matrix<T, 2, 1>& rM)
{
    using namespace boost::numeric::ublas;

    c_matrix<T, 1, 2> inverse;
    T det = Determinant(rM);

    inverse(0,0) = rM(0,0)/det/det;
    inverse(0,1) = rM(1,0)/det/det;

    return inverse;
}

/**
 * Calculates the pseudo-inverse of a 3x1 matrix.
 *
 * The pseudo inverse is given by pinv(T) = (T'T)^(-1)*T'
 */
template<class T>
boost::numeric::ublas::c_matrix<T, 1, 3> Inverse(const boost::numeric::ublas::c_matrix<T, 3, 1>& rM)
{
    using namespace boost::numeric::ublas;

    c_matrix<T, 1, 3> inverse;
    T det = Determinant(rM);

    inverse(0,0) = rM(0,0)/det/det;
    inverse(0,1) = rM(1,0)/det/det;
    inverse(0,2) = rM(2,0)/det/det;

    return inverse;
}

template<class T>
c_vector<T, 3> VectorProduct(const c_vector<T, 3>& rA, const c_vector<T, 3>& rB)
{
    /// This is a cross-product, only implemented for 3-vectors
    /// What do you get when you cross an elephant with a banana?
    /// \todo remove Joe's awful joke

    c_vector<T, 3> result;

    double x1 = rA(0);
    double y1 = rA(1);
    double z1 = rA(2);
    double x2 = rB(0);
    double y2 = rB(1);
    double z2 = rB(2);

    result(0) = y1*z2 - z1*y2;
    result(1) = z1*x2 - x1*z2;
    result(2) = x1*y2 - y1*x2;

    return result;
}

template<class T>
T Trace(const c_matrix<T, 1, 1>& rM)
{
    return rM(0,0);
}

template<class T>
T Trace(const c_matrix<T, 2, 2>& rM)
{
    return rM(0,0) + rM(1,1);
}

template<class T>
T Trace(const c_matrix<T, 3, 3>& rM)
{
    return rM(0,0) + rM(1,1) + rM(2,2);
}

template<class T>
T Trace(const c_matrix<T, 4, 4>& rM)
{
    return rM(0,0) + rM(1,1) + rM(2,2) + rM(3,3);
}

/**
 *  Second invariant of a matrix. Implementation only correct for
 *  a SYMMETRIC matrix though. It is up to the user to check the
 *  input matrix is symmetric.
 */
template<class T>
T SecondInvariant(const c_matrix<T, 3, 3>& rM)
{
    return    rM(0,0)*rM(1,1) + rM(1,1)*rM(2,2) + rM(2,2)*rM(0,0)
            - rM(1,0)*rM(1,0) - rM(2,1)*rM(2,1) - rM(2,0)*rM(2,0);
}

/**
 *  Second invariant of a 2d matrix, ie the determinant. This function
 *  is mainly here just so that the same code can be used in 2d and 3d.
 */
template<class T>
T SecondInvariant(const c_matrix<T, 2, 2>& rM)
{
    return Determinant(rM);
}

/**
 * Convenience functions for quickly creating test vectors.
 */
c_vector<double, 1> Create_c_vector(double x);

c_vector<double, 2> Create_c_vector(double x, double y);

c_vector<double, 3> Create_c_vector(double x, double y, double z);

/** Use LAPACK functionality to find the eigenvector corresponding
 *  real eigenvalue which is smallest in magnitude.
 * Caveat: if there are zero eigenvalues they are ignored.
 * It's the smallest magnitude non-zero real eigenvalue which is used.
 *
 * @param 3x3 matrix is question
 * @return 3-vector corresponding to right-eigenvector in question
 */
c_vector<double,3> CalculateEigenvectorForSmallestNonzeroEigenvalue(c_matrix<double, 3, 3>& rA);

/**
 * Replacement "pow" function
 * @param x
 * @param a  exponent
 * @return x^a, x**a.
 */
double SmallPow(double x, unsigned exponent);

/**
 * Uses fmod to determine if smallerNumber divides the largerNumber.
 * We expect smallerNumber/largerNumber <= 1 and therefore
 * fmod(largerNumber,smallerNumber) should be close zero or close to  smallerNumber
 * @param smallerNumber the smaller
 * @param largerNumber the larger
 */
bool Divides(double smallerNumber, double largerNumber);

#endif /*UBLASCUSTOMFUNCTIONS_HPP_*/

