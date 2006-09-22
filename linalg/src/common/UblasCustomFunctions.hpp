#ifndef UBLASCUSTOMFUNCTIONS_HPP_
#define UBLASCUSTOMFUNCTIONS_HPP_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>

#include <iostream>
#include <fstream>
#include <vector>

#define TINY 1E-13

using namespace boost::numeric::ublas;
/**
 * Get the determinant of a ublas matrix
 * N.B. Do not call with non-square matrix
 *
 *
 * @param m The matrix of which to find the determinant. Must be square.
 */
template<class T>
T Determinant(const boost::numeric::ublas::c_matrix<T,1,1> &m)
{
    using namespace boost::numeric::ublas;
    
    return m(0,0);
};
/**
 * Return the determinant of a submatrix after removing a particular row and column
 */
template<class T>
T SubDeterminant(const boost::numeric::ublas::c_matrix<T,1,1> &m, const int missrow, const int misscol)
{
    using namespace boost::numeric::ublas;
    
    assert(missrow==0);
    assert(misscol==0);
    return 1.0;
};
template<class T>
T Determinant(const boost::numeric::ublas::c_matrix<T,2,2> &m)
{
    using namespace boost::numeric::ublas;
    
    return m(0,0)*m(1,1)
           - m(1,0)*m(0,1);
};

/**
 * Return the determinant of a submatrix after removing a particular row and column
 */
template<class T>
T SubDeterminant(const boost::numeric::ublas::c_matrix<T,2,2> &m, const int missrow, const int misscol)
{
    using namespace boost::numeric::ublas;
    
    assert(missrow>=0 && missrow<2);
    assert(misscol>=0 && misscol<2);
    
    int row=(missrow==1)?0:1;
    int col=(misscol==1)?0:1;
    return m(row,col);
};

template<class T>
T Determinant(const boost::numeric::ublas::c_matrix<T,3,3> &m)
{
    using namespace boost::numeric::ublas;
    
    return   m(0,0) *
             (m(1,1)*m(2,2) - m(1,2)*m(2,1))
             - m(0,1) *
             (m(1,0)*m(2,2) - m(1,2)*m(2,0))
             + m(0,2) *
             (m(1,0)*m(2,1) - m(1,1)*m(2,0));
};
/**
 * Return the determinant of a submatrix after removing a particular row and column
 */
template<class T>
T SubDeterminant(const boost::numeric::ublas::c_matrix<T,3,3> &m, const int missrow, const int misscol)
{
    using namespace boost::numeric::ublas;
    
    assert(missrow>=0 && missrow<3);
    assert(misscol>=0 && misscol<3);
    
    int lorow=(missrow==0)?1:0;
    int hirow=(missrow==2)?1:2;
    int locol=(misscol==0)?1:0;
    int hicol=(misscol==2)?1:2;
    return (m(lorow,locol)*m(hirow,hicol) - m(lorow,hicol)*m(hirow,locol));
};

/**
 * Get the inverse of a ublas matrix
 * N.B. Do not call with non-square matrix
 *
 * @param m The matrix of which to find the inverse. Must be square.
 */
template<class T>
boost::numeric::ublas::c_matrix<T, 1, 1> Inverse(const boost::numeric::ublas::c_matrix<T, 1, 1> &m)
{
    using namespace boost::numeric::ublas;
    
    c_matrix<T, 1, 1> inverse;
    T det = Determinant(m);
    assert( fabs(det) > TINY ); //Else it is a singular matrix
    inverse(0,0) =  1.0/det;
    return inverse;
};

template<class T>
boost::numeric::ublas::c_matrix<T, 2, 2> Inverse(const boost::numeric::ublas::c_matrix<T, 2, 2> &m)
{
    using namespace boost::numeric::ublas;
    
    c_matrix<T, 2, 2> inverse;
    T det = Determinant(m);
    assert( fabs(det) > TINY ); //Else it is a singular matrix
    inverse(0,0)  =  m(1,1)/det;
    inverse(0,1)  = -m(0,1)/det;
    inverse(1,0)  = -m(1,0)/det;
    inverse(1,1)  =  m(0,0)/det;
    return inverse;
};

template<class T>
boost::numeric::ublas::c_matrix<T, 3, 3> Inverse(const boost::numeric::ublas::c_matrix<T, 3, 3> &m)
{
    using namespace boost::numeric::ublas;
    
    c_matrix<T, 3, 3> inverse;
    T det = Determinant(m);
    assert( fabs(det) > TINY ); //Else it is a singular matrix
    inverse(0,0)  =  (m(1,1)*m(2,2)-m(1,2)*m(2,1))/det;
    inverse(1,0)  =  -(m(1,0)*m(2,2)-m(1,2)*m(2,0))/det;
    inverse(2,0)  =  (m(1,0)*m(2,1)-m(1,1)*m(2,0))/det;
    inverse(0,1)  =  -(m(0,1)*m(2,2)-m(0,2)*m(2,1))/det;
    inverse(1,1)  =  (m(0,0)*m(2,2)-m(0,2)*m(2,0))/det;
    inverse(2,1)  =  -(m(0,0)*m(2,1)-m(0,1)*m(2,0))/det;
    inverse(0,2)  =  (m(0,1)*m(1,2)-m(0,2)*m(1,1))/det;
    inverse(1,2)  = - (m(0,0)*m(1,2)-m(0,2)*m(1,0))/det;
    inverse(2,2)  =  (m(0,0)*m(1,1)-m(0,1)*m(1,0))/det;
    return inverse;
};

template<class T>
c_vector<T, 3> VectorProduct(const c_vector<T, 3> &a, const c_vector<T, 3> &b)
{
    //This is a cross-product
    // What do you get when you cross an elephant with a banana?
    //only implemented for 3-vectors
    
    c_vector<T, 3>  result;
    
    double x1=a(0);
    double y1=a(1);
    double z1=a(2);
    double x2=b(0);
    double y2=b(1);
    double z2=b(2);
    
    result(0)=y1*z2-z1*y2;
    result(1)=z1*x2-x1*z2;
    result(2)=x1*y2-y1*x2;
    
    return result;
}

#endif /*UBLASCUSTOMFUNCTIONS_HPP_*/
