#ifndef UBLASCUSTOMFUNCTIONS_HPP_
#define UBLASCUSTOMFUNCTIONS_HPP_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>

#include <iostream>
#include <fstream>
#include <vector>

#define TINY 1E-13

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

template<class T>
T Determinant(const boost::numeric::ublas::c_matrix<T,2,2> &m)
{
    using namespace boost::numeric::ublas;
    
    return m(0,0)*m(1,1) 
           - m(1,0)*m(0,1);
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

#endif /*UBLASCUSTOMFUNCTIONS_HPP_*/
