#ifndef UBLASCUSTOMFUNCTIONS_HPP_
#define UBLASCUSTOMFUNCTIONS_HPP_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>

#include <iostream>
#include <fstream>
#include <vector>

#define TINY 0.00000001

/**
 * Get the determinant of a ublas matrix
 * N.B. Do not call with non-square matrix
 * 
 * 
 * @param m The matrix of which to find the determinant. Must be square.
 */
template<class T, unsigned rows, unsigned cols>
T Determinant(const boost::numeric::ublas::c_matrix<T,rows,cols> &m)
{
    using namespace boost::numeric::ublas;
    
    
    // Matrix must be square and of size 1, 2 or 3
    assert(rows == cols);
    assert(rows>0 && rows<4);
     
    T det;
    switch( rows )
    {
        case 1 :
            det = m(0,0);
            break;
        case 2 :
            det = m(0,0)*m(1,1) 
                - m(1,0)*m(0,1);
            break;
        case 3 :
            det = m(0,0) * 
            (m(1,1)*m(2,2) - m(1,2)*m(2,1))
             - m(0,1) * 
             (m(1,0)*m(2,2) - m(1,2)*m(2,0))
             + m(0,2) * 
             (m(1,0)*m(2,1) - m(1,1)*m(2,0));
    }
    return det;
};
/**
 *
 * N.B. Do not call with non-square matrix
 * 
 * @param m The matrix of which to find the inverse. Must be square.
 */

template<class T, unsigned rows, unsigned cols>
boost::numeric::ublas::c_matrix<T, rows, cols> Inverse(const boost::numeric::ublas::c_matrix<T,rows,cols> &m)
{
    using namespace boost::numeric::ublas;
        
    // Matrix must be square and of size 1, 2 or 3
    assert(rows == cols);
    assert(rows>0 && rows<4);
    c_matrix<T, rows, cols> inverse;
    T det = Determinant(m);
    assert( fabs(det) > TINY ); //Else it is a singular matrix
    
    switch( rows )
    {
        case 1 :
            inverse(0,0) =  1.0/det;
            break;
        case 2 :
            inverse(0,0)  =  m(1,1)/det;
            inverse(0,1)  = -m(0,1)/det;
            inverse(1,0)  = -m(1,0)/det;
            inverse(1,1)  =  m(0,0)/det;
            break;
        case 3 :
            inverse(0,0)  =  (m(1,1)*m(2,2)-m(1,2)*m(2,1))/det;
            inverse(1,0)  =  -(m(1,0)*m(2,2)-m(1,2)*m(2,0))/det;
            inverse(2,0)  =  (m(1,0)*m(2,1)-m(1,1)*m(2,0))/det;
            inverse(0,1)  =  -(m(0,1)*m(2,2)-m(0,2)*m(2,1))/det;
            inverse(1,1)  =  (m(0,0)*m(2,2)-m(0,2)*m(2,0))/det;
            inverse(2,1)  =  -(m(0,0)*m(2,1)-m(0,1)*m(2,0))/det;
            inverse(0,2)  =  (m(0,1)*m(1,2)-m(0,2)*m(1,1))/det;
            inverse(1,2)  = - (m(0,0)*m(1,2)-m(0,2)*m(1,0))/det;
            inverse(2,2)  =  (m(0,0)*m(1,1)-m(0,1)*m(1,0))/det;
    }
    return inverse;   
         
    
};

#endif /*UBLASCUSTOMFUNCTIONS_HPP_*/
