#ifndef _VECTORDOUBLE_HPP_
#define _VECTORDOUBLE_HPP_

#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/io.hpp>
#include <cassert>
#include <iostream>
#include <math.h>

using namespace boost::numeric::ublas;

class VectorDouble
{
    private:
       int mSize;
       
       c_vector<double,1> *mpVectorOf1;
       c_vector<double,2> *mpVectorOf2;
       c_vector<double,3> *mpVectorOf3;
       c_vector<double,4> *mpVectorOf4;
      
    public:
        //VectorDouble(); // For use in alocating a Double vector 
        VectorDouble(int Size);
        ~VectorDouble();

          VectorDouble(const VectorDouble& rOtherVector);
          void operator=(const VectorDouble& rOtherVector);
          VectorDouble operator+(const VectorDouble& rSomeVector1);
          VectorDouble operator-(const VectorDouble& rSomeVector1);
          VectorDouble operator*(double Scalar);

          double& VectorDouble::operator()(int Entry) const;

          int Size( void ) const;
          double dot(const VectorDouble& rOtherVector) const;
          void ResetToZero( void );
          VectorDouble VectorProduct(const VectorDouble& rOtherVector);
          friend VectorDouble operator*(double Scalar, const VectorDouble& rSomeVector);
          double L2Norm( void );
          
          c_vector<double, 1>* GetUblasHandle1( void ) const;
          c_vector<double, 2>* GetUblasHandle2( void ) const;
          c_vector<double, 3>* GetUblasHandle3( void ) const;
          c_vector<double, 4>* GetUblasHandle4( void ) const;
          
          c_vector<double, 1>& rGetUblasHandle1( void );
          c_vector<double, 2>& rGetUblasHandle2( void );
          c_vector<double, 3>& rGetUblasHandle3( void );
          c_vector<double, 4>& rGetUblasHandle4( void );
};



#endif //_VECTORDOUBLE_HPP_
