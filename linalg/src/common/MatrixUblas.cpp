#include "MatrixUblas.hpp"
#include "Exception.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <cassert>
#include <iostream>
#include <math.h>

#define TINY 0.00000001

MatrixUblas::MatrixUblas(int numRows, int numColumns)
{
    assert(numRows == numColumns );
    
    mSize = numRows;
    
    assert(mSize > 0);
    
    using namespace boost::numeric::ublas;
    
    switch (mSize)
    {
        case 1:
            mpMatrixOf1 = new c_matrix<double,1,1>(1,1); 
            break;
        case 2:
            mpMatrixOf2 = new c_matrix<double,2,2>(2,2); 
            break; 
        case 3:
            mpMatrixOf3 = new c_matrix<double,3,3>(3,3); 
            break; 
        case 4:
            mpMatrixOf4 = new c_matrix<double,4,4>(4,4); 
            break; 
        
        default:
            // Matrix biger than  4 Throw Exception
            throw("Matrix size larger than 4");
            break;
    }
    
     
}

MatrixUblas::MatrixUblas(const MatrixUblas& rOtherMatrix)
{
    mSize = rOtherMatrix.mSize;
    
    switch (mSize)
    {
        case 1:
             mpMatrixOf1 = new c_matrix<double,1,1>(*(rOtherMatrix.mpMatrixOf1)); 
             break;
        case 2:
             mpMatrixOf2 = new c_matrix<double,2,2>(*(rOtherMatrix.mpMatrixOf2)); 
             break; 
        case 3:
             mpMatrixOf3 = new c_matrix<double,3,3>(*(rOtherMatrix.mpMatrixOf3)); 
             break; 
        case 4:
             mpMatrixOf4 = new c_matrix<double,4,4>(*(rOtherMatrix.mpMatrixOf4)); 
             break; 
        
        default:
            // Matrix biger than  4 Throw Exception
            throw("Matrix size larger than 4");
            break;
    }
    
    
}


MatrixUblas::~MatrixUblas()
{
       switch (mSize)
    {
        case 1:
            delete (mpMatrixOf1);
            break;
        case 2:
            delete (mpMatrixOf2);
            break;
        case 3:
            delete (mpMatrixOf3);
            break;
        case 4:
            delete (mpMatrixOf4);
            break;
        default:
            // Matrix biger than  4 Throw Exception
            throw("Matrix size larger than 4");
            break;
    }     
}



double& MatrixUblas::operator()(int Row, int Column) const
{
    assert(Row > -1);
    assert(Row < mSize);
    assert(Column > -1);
    assert(Column < mSize);
    
      switch(mSize) {
        case 1:
            return (*mpMatrixOf1) (Row,Column);
        case 2:
            return (*mpMatrixOf2) (Row,Column);
        case 3:
            return (*mpMatrixOf3) (Row,Column);
        case 4:
            return (*mpMatrixOf4) (Row,Column);
        default:
            // Matrix biger than  4 Throw Exception
            throw("Matrix size larger than 4");
    }

}

MatrixUblas& MatrixUblas::operator=(const MatrixUblas& rOtherMatrix)
{
    assert( mSize == rOtherMatrix.mSize);
    
    switch(mSize) {
        case 1:
            *mpMatrixOf1=*(rOtherMatrix.mpMatrixOf1);
            break;
        case 2:
            *mpMatrixOf2=*(rOtherMatrix.mpMatrixOf2);
            break;
        case 3:
            *mpMatrixOf3=*(rOtherMatrix.mpMatrixOf3);
            break;
        case 4:
            *mpMatrixOf4=*(rOtherMatrix.mpMatrixOf4);
            break;
      default:
            // Matrix biger than  4 Throw Exception
            throw("Matrix size larger than 4");
    }
           
}

MatrixUblas MatrixUblas::operator*(double scalar)
{
    MatrixUblas result(mSize, mSize);
    for (int i=0; i<mSize; i++)
    {
        for (int j=0; j<mSize; j++)
        { 
           result(i,j) = scalar * (*this)(i,j);
        }
    }
    return result;
}

MatrixUblas MatrixUblas::Identity(int Size)
{
    MatrixUblas Eye(Size, Size);
    for (int i=0; i<Size; i++)
    {
        Eye(i,i) = 1.0;
    }
    return Eye;
}
 
 int MatrixUblas::Rows( void ) const
 {
    return mSize;
 }
 
 int MatrixUblas::Columns( void ) const
 {
    return mSize;
 }
 
 double MatrixUblas::Determinant( void ) const
 {
    assert( mSize > 0 && mSize < 4);
    double det = 0.0;
    switch( mSize )
    {
        case 1 :
            det = (*mpMatrixOf1)(0,0);
            break;
        case 2 :
            det = (*mpMatrixOf2)(0,0)*(*mpMatrixOf2)(1,1) 
                - (*mpMatrixOf2)(1,0)*(*mpMatrixOf2)(0,1);
            break;
        case 3 :
            det = (*mpMatrixOf3)(0,0) * 
            ((*mpMatrixOf3)(1,1)*(*mpMatrixOf3)(2,2) - (*mpMatrixOf3)(1,2)*(*mpMatrixOf3)(2,1))
             - (*mpMatrixOf3)(0,1) * 
             ((*mpMatrixOf3)(1,0)*(*mpMatrixOf3)(2,2) - (*mpMatrixOf3)(1,2)*(*mpMatrixOf3)(2,0))
             + (*mpMatrixOf3)(0,2) * 
             ((*mpMatrixOf3)(1,0)*(*mpMatrixOf3)(2,1) - (*mpMatrixOf3)(1,1)*(*mpMatrixOf3)(2,0));
    }
    return det;
    // If this were to be implemented for 4x4 or larger matrices then
    // the correct thing to do is to make an LU factorisation (Ublas will
    // do this for us), and then multiply the diagonal elements
}
 

 MatrixUblas MatrixUblas::Inverse( void ) const
 {
    assert( mSize > 0 && mSize < 4);
    MatrixUblas Inverse(mSize,mSize);
    double Det = Determinant();
    assert( fabs(Det) > TINY ); //Else it is a singular matrix
    switch( mSize )
    {
        case 1 :
            Inverse(0,0) =  1.0/Det;
            break;
        case 2 :
            Inverse(0,0)  =  (*mpMatrixOf2)(1,1)/Det;
            Inverse(0,1)  = -(*mpMatrixOf2)(0,1)/Det;
            Inverse(1,0)  = -(*mpMatrixOf2)(1,0)/Det;
            Inverse(1,1)  =  (*mpMatrixOf2)(0,0)/Det;
            break;
        case 3 :
            Inverse(0,0)  =  ((*mpMatrixOf3)(1,1)*(*mpMatrixOf3)(2,2)-(*mpMatrixOf3)(1,2)*(*mpMatrixOf3)(2,1))/Det;
            Inverse(1,0)  =  -((*mpMatrixOf3)(1,0)*(*mpMatrixOf3)(2,2)-(*mpMatrixOf3)(1,2)*(*mpMatrixOf3)(2,0))/Det;
            Inverse(2,0)  =  ((*mpMatrixOf3)(1,0)*(*mpMatrixOf3)(2,1)-(*mpMatrixOf3)(1,1)*(*mpMatrixOf3)(2,0))/Det;
            Inverse(0,1)  =  -((*mpMatrixOf3)(0,1)*(*mpMatrixOf3)(2,2)-(*mpMatrixOf3)(0,2)*(*mpMatrixOf3)(2,1))/Det;
            Inverse(1,1)  =  ((*mpMatrixOf3)(0,0)*(*mpMatrixOf3)(2,2)-(*mpMatrixOf3)(0,2)*(*mpMatrixOf3)(2,0))/Det;
            Inverse(2,1)  =  -((*mpMatrixOf3)(0,0)*(*mpMatrixOf3)(2,1)-(*mpMatrixOf3)(0,1)*(*mpMatrixOf3)(2,0))/Det;
            Inverse(0,2)  =  ((*mpMatrixOf3)(0,1)*(*mpMatrixOf3)(1,2)-(*mpMatrixOf3)(0,2)*(*mpMatrixOf3)(1,1))/Det;
            Inverse(1,2)  = - ((*mpMatrixOf3)(0,0)*(*mpMatrixOf3)(1,2)-(*mpMatrixOf3)(0,2)*(*mpMatrixOf3)(1,0))/Det;
            Inverse(2,2)  =  ((*mpMatrixOf3)(0,0)*(*mpMatrixOf3)(1,1)-(*mpMatrixOf3)(0,1)*(*mpMatrixOf3)(1,0))/Det;
    }
    return Inverse;   
    
 }
 
 
 
 VectorUblas MatrixUblas::operator*(const VectorUblas& rSomeVector) const
{
    assert(mSize==rSomeVector.Size());
    VectorUblas result(mSize);
    
   
    
    for( int i = 0; i < mSize; i++)
    {
        for(int j = 0; j < mSize; j++)
        {
            result(i) += (*this)(i,j) * rSomeVector(j);
        }
    }
    return result;
}
 
 
 
 
