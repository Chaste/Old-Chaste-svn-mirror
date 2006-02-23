#ifndef _VECTORUBLAS_HPP_
#define _VECTORUBLAS_HPP_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <cassert>
#include <iostream>
#include <math.h>

using namespace boost::numeric::ublas;

class VectorUblas
{
    private:
       int mSize;
       
       c_vector<double,1> *mpVectorOf1;
       c_vector<double,2> *mpVectorOf2;
       c_vector<double,3> *mpVectorOf3;
       c_vector<double,4> *mpVectorOf4;
      
    public:
        //VectorUblas(); // For use in alocating a Ublas vector 
        VectorUblas(int Size);
        ~VectorUblas();

          VectorUblas(const VectorUblas& rOtherVector);
          VectorUblas& operator=(const VectorUblas& rOtherVector);
          VectorUblas operator+(const VectorUblas& rSomeVector1);
          VectorUblas operator-(const VectorUblas& rSomeVector1);
          VectorUblas operator*(double Scalar);

          double& VectorUblas::operator()(int Entry) const;

          int Size( void ) const;
          double dot(const VectorUblas& rOtherVector) const;
          void ResetToZero( void );
          VectorUblas VectorProduct(const VectorUblas& rOtherVector);
          friend VectorUblas operator*(double Scalar, const VectorUblas& rSomeVector);
          double L2Norm( void );

};



#endif //_VECTORUBLAS_HPP_
