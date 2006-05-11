#ifndef MATRIXDOUBLE_HPP_
#define MATRIXDOUBLE_HPP_

#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/io.hpp>
#include <cassert>
#include <iostream>
#include <math.h>
#include "VectorDouble.hpp"

using namespace boost::numeric::ublas;

class MatrixDouble
{
    private:
        int mSize; //Only square matrices =>  mRows == mColumns
     
        /**
         * MatrixDouble is implemented using ublas matrices.
         */
        c_matrix<double,1,1> *mpMatrixOf1;
        c_matrix<double,2,2> *mpMatrixOf2;
        c_matrix<double,3,3> *mpMatrixOf3;
        c_matrix<double,4,4> *mpMatrixOf4;
        
        int mNumberOfElements;
                
    public:
          MatrixDouble(int Rows, int Columns);
   
          MatrixDouble(const MatrixDouble& rOtherMatrix);
          ~MatrixDouble();
        
          void operator=(const MatrixDouble& rOtherMatrix);
          double &operator()(int Row, int Column) const;
          static MatrixDouble Identity(int Size);
          int Rows( void ) const;
          int Columns( void ) const;
          MatrixDouble Inverse( void ) const;
          double Determinant( void ) const;
          void ResetToZero( void );
//        
          MatrixDouble operator*(double scalar);
          friend MatrixDouble operator*(const double scalar, const MatrixDouble& rMatrix);
          friend MatrixDouble operator*(const MatrixDouble &rLeftMatrix, const MatrixDouble &rRightMatrix);
          friend MatrixDouble operator+(const MatrixDouble &rLeftMatrix, const MatrixDouble &rRightMatrix);
          friend MatrixDouble operator-(const MatrixDouble &rLeftMatrix, const MatrixDouble &rRightMatrix);
//        
//        void VectorPostMultiply(const VectorDouble& rOperandVector, VectorDouble& rResultVector) const;
//        
//        
          friend VectorDouble operator*(const VectorDouble& rSomeVector, const MatrixDouble& rSomeMatrix);
          VectorDouble operator*(const VectorDouble& rSomeVector) const;
          MatrixDouble Transpose() const;
          bool IsSquare( void ) const;
          double GetTrace() const;
          double GetFirstInvariant() const;
          double GetSecondInvariant() const;
          double GetThirdInvariant() const;

          c_matrix<double, 1,1>* GetUblasHandle1( void ) const;
          c_matrix<double, 2,2>* GetUblasHandle2( void ) const;
          c_matrix<double, 3,3>* GetUblasHandle3( void ) const;
          c_matrix<double, 4,4>* GetUblasHandle4( void ) const;

          c_matrix<double, 1,1>& rGetUblasHandle1( void );
          c_matrix<double, 2,2>& rGetUblasHandle2( void );
          c_matrix<double, 3,3>& rGetUblasHandle3( void );
          c_matrix<double, 4,4>& rGetUblasHandle4( void );

};

#endif /*MATRIXDOUBLE_HPP_*/
