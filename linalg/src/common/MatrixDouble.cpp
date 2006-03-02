#include "MatrixDouble.hpp"
#include "Exception.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <cassert>
#include <iostream>
#include <math.h>

#define TINY 0.00000001

MatrixDouble::MatrixDouble(int numRows, int numColumns)
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

MatrixDouble::MatrixDouble(const MatrixDouble& rOtherMatrix)
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


MatrixDouble::~MatrixDouble()
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



double& MatrixDouble::operator()(int Row, int Column) const
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

MatrixDouble& MatrixDouble::operator=(const MatrixDouble& rOtherMatrix)
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

MatrixDouble MatrixDouble::operator*(double scalar)
{
    MatrixDouble result(mSize, mSize);
    for (int i=0; i<mSize; i++)
    {
        for (int j=0; j<mSize; j++)
        { 
           result(i,j) = scalar * (*this)(i,j);
        }
    }
    return result;
}

MatrixDouble MatrixDouble::Identity(int Size)
{
    MatrixDouble Eye(Size, Size);
    for (int i=0; i<Size; i++)
    {
        Eye(i,i) = 1.0;
    }
    return Eye;
}
 
 int MatrixDouble::Rows( void ) const
 {
    return mSize;
 }
 
 int MatrixDouble::Columns( void ) const
 {
    return mSize;
 }
 
 double MatrixDouble::Determinant( void ) const
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
    // the correct thing to do is to make an LU factorisation (Double will
    // do this for us), and then multiply the diagonal elements
}
 

 MatrixDouble MatrixDouble::Inverse( void ) const
 {
    assert( mSize > 0 && mSize < 4);
    MatrixDouble Inverse(mSize,mSize);
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
 
 
 
 VectorDouble MatrixDouble::operator*(const VectorDouble& rSomeVector) const
{
    assert(mSize==rSomeVector.Size());
    VectorDouble result(mSize);
    
   
    
    for( int i = 0; i < mSize; i++)
    {
        for(int j = 0; j < mSize; j++)
        {
            result(i) += (*this)(i,j) * rSomeVector(j);
        }
    }
    return result;
}
 
 MatrixDouble MatrixDouble::Transpose() const
{
    MatrixDouble result(mSize, mSize);
    for( int i = 0; i < mSize; i++)
    {
        for(int j = 0; j < mSize; j++)
        {
            result(j,i) = (*this)(i,j);
        }
    }
    return result;
}

VectorDouble operator* (const VectorDouble& rSomeVector, const MatrixDouble& rSomeMatrix)
{
    assert(rSomeVector.Size()==rSomeMatrix.Rows());
    VectorDouble result(rSomeMatrix.Columns());
    int index;
    for( int i = 0; i < rSomeMatrix.Columns(); i++)
    {
        for(int j = 0; j < rSomeMatrix.Rows(); j++)
        {
            result(i) += rSomeMatrix(j,i)*rSomeVector(j);
        }
    }
    return result;
}
 
 void MatrixDouble::ResetToZero( void )
{
    for( int i = 0; i < mSize; i++)
        {
            for( int j = 0; j < mSize; j++)
            {
                (*this)(i,j) = 0.0;
            }
        }
}

double MatrixDouble::GetTrace() const
{
    double sum=0;
    for(int i=0; i<mSize; i++)
    {
        sum += (*this)(i,i);
    }
    return sum;
}


double MatrixDouble::GetFirstInvariant() const
{
    assert(mSize<=3);
    
    return GetTrace();
}

double MatrixDouble::GetSecondInvariant() const
{
    assert(mSize<=3);
    assert(mSize>1);
    
    double ret;
    if(mSize==2)
    {
        // second invariant in 2d is the determinant
        ret = Determinant();
    }
    else if(mSize==3)
    {
        // second invariant in 3d is 0.5*(tr(C^2)-tr(C)^2;
        // ie C[0][1]*C[1][0] + C[1][2]*C[2][1] + C[2][0]*C[0][2] - C[0][0]*C[1][1] - C[1][1]*C[2][2] - C[2][2]*C[0][0];
        
        ret =   (*mpMatrixOf3)(0,1)*(*mpMatrixOf3)(1,0) + (*mpMatrixOf3)(0,2)*(*mpMatrixOf3)(2,0) 
              + (*mpMatrixOf3)(1,2)*(*mpMatrixOf3)(2,1) - (*mpMatrixOf3)(0,0)*(*mpMatrixOf3)(1,1)
              - (*mpMatrixOf3)(1,1)*(*mpMatrixOf3)(2,2) - (*mpMatrixOf3)(2,2)*(*mpMatrixOf3)(0,0);
    }
    else
    {
        // shouldn't be possible to get here
        assert(0);
    }
    
    return ret;
}
    
double MatrixDouble::GetThirdInvariant() const
{
    assert(mSize==3);

    return Determinant();
}   
 
 /// \todo make this more efficient?
 // We should be using a Double function, if there is one available.
MatrixDouble operator*(const MatrixDouble &rLeftMatrix, const MatrixDouble &rRightMatrix)
{
    assert(rLeftMatrix.Columns()==rRightMatrix.Rows());

    MatrixDouble result(rLeftMatrix.Rows(), rRightMatrix.Columns());
    
    for(int i=0; i<rLeftMatrix.Rows(); i++)
    {
        for(int j=0; j<rRightMatrix.Columns(); j++)
        {
            for(int k=0; k<rLeftMatrix.Columns(); k++)
            {
                result(i,j) += rLeftMatrix(i,k)*rRightMatrix(k,j);
            }
        }
    }
    return result;
}

MatrixDouble operator+(const MatrixDouble &rLeftMatrix, const MatrixDouble &rRightMatrix)
{
    assert(rLeftMatrix.Rows()==rRightMatrix.Rows());
    assert(rLeftMatrix.Columns()==rRightMatrix.Columns());

    MatrixDouble result(rLeftMatrix.Rows(), rLeftMatrix.Columns());
    
    for(int i=0; i<rLeftMatrix.Rows(); i++)
    {
        for(int j=0; j<rLeftMatrix.Columns(); j++)
        {
            result(i,j) += rLeftMatrix(i,j) + rRightMatrix(i,j);
        }
    }
    return result;
}


/// \todo make this more efficient?
MatrixDouble operator-(const MatrixDouble &rLeftMatrix, const MatrixDouble &rRightMatrix)
{
    assert(rLeftMatrix.Rows()==rRightMatrix.Rows());
    assert(rLeftMatrix.Columns()==rRightMatrix.Columns());

    MatrixDouble result(rLeftMatrix.Rows(), rLeftMatrix.Columns());
    
    for(int i=0; i<rLeftMatrix.Rows(); i++)
    {
        for(int j=0; j<rLeftMatrix.Columns(); j++)
        {
            result(i,j) += rLeftMatrix(i,j) - rRightMatrix(i,j);
        }
    }
    return result;
}

MatrixDouble operator*(double scalar, const MatrixDouble &rMatrix)
{
    MatrixDouble result(rMatrix.Rows(), rMatrix.Columns());
    for (int i=0; i<rMatrix.mSize; i++)
    {
        for (int j=0; j<rMatrix.mSize; j++)
        { 
           result(i,j) = scalar * rMatrix(i,j);
        }
    }
    return result;
}

bool MatrixDouble::IsSquare() const
{
    return true;
}




