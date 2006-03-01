#ifndef MATRIXUBLAS_HPP_
#define MATRIXUBLAS_HPP_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <cassert>
#include <iostream>
#include <math.h>
#include "VectorUblas.hpp"

using namespace boost::numeric::ublas;

class MatrixUblas
{
    private:
        int mSize; //Only square matrices =>  mRows == mColumns
     
        c_matrix<double,1,1> *mpMatrixOf1;
        c_matrix<double,2,2> *mpMatrixOf2;
        c_matrix<double,3,3> *mpMatrixOf3;
        c_matrix<double,4,4> *mpMatrixOf4;
        int mNumberOfElements;
                
    public:
          MatrixUblas(int Rows, int Columns);
   
          MatrixUblas(const MatrixUblas& rOtherMatrix);
          ~MatrixUblas();
        
          MatrixUblas& operator=(const MatrixUblas& rOtherMatrix);
          double &operator()(int Row, int Column) const;
          static MatrixUblas Identity(int Size);
          int Rows( void ) const;
          int Columns( void ) const;
          MatrixUblas Inverse( void ) const;
          double Determinant( void ) const;
          void ResetToZero( void );
//        
          MatrixUblas operator*(double scalar);
          friend MatrixUblas operator*(const double scalar, const MatrixUblas& rMatrix);
          friend MatrixUblas operator*(const MatrixUblas &rLeftMatrix, const MatrixUblas &rRightMatrix);
          friend MatrixUblas operator+(const MatrixUblas &rLeftMatrix, const MatrixUblas &rRightMatrix);
          friend MatrixUblas operator-(const MatrixUblas &rLeftMatrix, const MatrixUblas &rRightMatrix);
//        
//        void VectorPostMultiply(const VectorDouble& rOperandVector, VectorDouble& rResultVector) const;
//        
//        
          friend VectorUblas operator*(const VectorUblas& rSomeVector, const MatrixUblas& rSomeMatrix);
          VectorUblas operator*(const VectorUblas& rSomeVector) const;
          MatrixUblas Transpose() const;
          bool IsSquare( void ) const;
          double GetTrace() const;
          double GetFirstInvariant() const;
          double GetSecondInvariant() const;
          double GetThirdInvariant() const;


};

#endif /*MATRIXUBLAS_HPP_*/
