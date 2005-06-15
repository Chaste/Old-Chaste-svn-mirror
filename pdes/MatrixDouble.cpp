#include "MatrixDouble.hpp"
#include <cassert>
#define TINY 0.00000001
#include <cmath>

#include <iostream>

MatrixDouble::MatrixDouble(int Rows, int Columns)
{
	assert(Rows > 0);
	assert(Columns > 0);
	mRows = Rows;
	mColumns = Columns;
	mNumberOfElements = mRows * mColumns;
	mElementArray = new double [mNumberOfElements];
	
	for (int i=0; i<mNumberOfElements; i++)
	{
		mElementArray[i] = 0.0;
	}
}


MatrixDouble::MatrixDouble(const MatrixDouble& rOtherMatrix)
{
	mRows = rOtherMatrix.mRows;
	mColumns = rOtherMatrix.mColumns;
	mNumberOfElements = mRows * mColumns;
	mElementArray = new double [mNumberOfElements];
	
	for (int i=0; i<mNumberOfElements; i++)
	{
		mElementArray[i] = rOtherMatrix.mElementArray[i];
	}
}


MatrixDouble::~MatrixDouble()
{
	delete[] mElementArray;
}


MatrixDouble& MatrixDouble::operator=(const MatrixDouble& rOtherMatrix)
{
	assert(mRows == rOtherMatrix.mRows);
	assert(mColumns == rOtherMatrix.mColumns);
	for (int i=0; i<mNumberOfElements; i++)
	{
		mElementArray[i] = rOtherMatrix.mElementArray[i];
	}
	return *this;
}


double& MatrixDouble::operator()(int Row, int Column) const
{
	assert(Row > -1);
	assert(Row < mRows);
	assert(Column > -1);
	assert(Column < mColumns);
	int Index = Column + mColumns*Row;
	return mElementArray[Index];
}


int MatrixDouble::Rows( void ) const
{
	return mRows;
}


int MatrixDouble::Columns( void ) const
{
	return mColumns;
}		


MatrixDouble& MatrixDouble::operator*(double scalar)
{
	for (int i=0; i<mNumberOfElements; i++)
	{
		mElementArray[i] *= scalar;
	}
	return *this;
}


MatrixDouble operator*(double scalar, const MatrixDouble &rMatrix) const
{
	MatrixDouble result(rMatrix.Rows(), rMatrix.Columns());
	for (int i=0; i<rMatrix.mNumberOfElements; i++)
	{
		result.mElementArray[i] = rMatrix.mElementArray[i] * scalar;
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


double MatrixDouble::Determinant() const
{
	assert( mRows == mColumns );
	assert( mRows > 0 && mRows < 4);
	double ret = 0.0;
	switch( mRows )
	{
		case 1 :
			ret = mElementArray[0];
			break;
		case 2 :
			ret = mElementArray[0]*mElementArray[3]-mElementArray[1]*mElementArray[2];
			break;
		case 3 :
			ret = mElementArray[0]*(mElementArray[8]*mElementArray[4] - mElementArray[5]*mElementArray[7])
             - mElementArray[3]*(mElementArray[1]*mElementArray[8] - mElementArray[2]*mElementArray[7])
             + mElementArray[6]*(mElementArray[1]*mElementArray[5] - mElementArray[2]*mElementArray[4]);
	}
	return ret;
}


MatrixDouble MatrixDouble::Inverse( void ) const
{
	assert( mRows == mColumns );
	assert( mRows > 0 && mRows < 4);
	MatrixDouble Inverse(mRows,mColumns);
	double Det = Determinant();
	assert( fabs(Det) > TINY );	
	switch( mRows )
	{
		case 1 :
			Inverse(0,0) =  1.0/Det;
			break;
		case 2 :
			Inverse(0,0)  =  mElementArray[3]/Det;
        	Inverse(0,1)  = -mElementArray[1]/Det;
        	Inverse(1,0)  = -mElementArray[2]/Det;
        	Inverse(1,1)  =  mElementArray[0]/Det;
        	break;
        case 3 :
        	Inverse(0,0)  =  (mElementArray[4]*mElementArray[8]-mElementArray[5]*mElementArray[7])/Det;
        	Inverse(1,0)  =  -(mElementArray[3]*mElementArray[8]-mElementArray[5]*mElementArray[6])/Det;
        	Inverse(2,0)  =  (mElementArray[3]*mElementArray[7]-mElementArray[4]*mElementArray[6])/Det;
        	Inverse(0,1)  =  -(mElementArray[1]*mElementArray[8]-mElementArray[2]*mElementArray[7])/Det;
        	Inverse(1,1)  =  (mElementArray[0]*mElementArray[8]-mElementArray[2]*mElementArray[6])/Det;
        	Inverse(2,1)  =  -(mElementArray[0]*mElementArray[7]-mElementArray[1]*mElementArray[6])/Det;
        	Inverse(0,2)  =  (mElementArray[1]*mElementArray[5]-mElementArray[2]*mElementArray[4])/Det;
        	Inverse(1,2)  = - (mElementArray[0]*mElementArray[5]-mElementArray[2]*mElementArray[3])/Det;
        	Inverse(2,2)  =  (mElementArray[0]*mElementArray[4]-mElementArray[1]*mElementArray[3])/Det;
	}
	return Inverse;
}


VectorDouble MatrixDouble::operator*(const VectorDouble& rSomeVector) const
{
	assert(mColumns==rSomeVector.Size());
	VectorDouble result(mRows);
	int index;
	for( int i = 0; i < mRows; i++)
	{
		for(int j = 0; j < mColumns; j++)
		{
			index = j + mColumns*i;
			result(i) += mElementArray[index]*rSomeVector(j);
		}
	}
	return result;
}

///
/// \todo make this more efficient?
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


///
/// \todo make this more efficient?
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

///
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

	
MatrixDouble MatrixDouble::Transpose() const
{
	MatrixDouble result(mColumns, mRows);
	for( int i = 0; i < mRows; i++)
	{
		for(int j = 0; j < mColumns; j++)
		{
			result(j,i) = (*this)(i,j);
		}
	}
	return result;
}

void MatrixDouble::ResetToZero( void )
{
	for (int i=0; i<mNumberOfElements; i++)
	{
		mElementArray[i]=0;
	}
}

bool MatrixDouble::IsSquare( void ) const
{
	return (mRows==mColumns);
}


double MatrixDouble::GetTrace() const
{
	assert(mRows==mColumns);
	double sum=0;
	for(int i=0; i<mRows; i++)
	{
		int index = i + mColumns*i;
		sum += mElementArray[index];
	}
	return sum;
}


double MatrixDouble::GetFirstInvariant() const
{
	assert(mRows==mColumns);
	assert(mRows<=3);
	
	return GetTrace();
}

double MatrixDouble::GetSecondInvariant() const
{
	assert(mRows==mColumns);
	assert(mRows<=3);
	assert(mRows>1);
	
	double ret;
	if(mRows==2)
	{
		// second invariant in 2d is the determinant
		ret = Determinant();
	}
	else if(mRows==3)
	{
		// second invariant in 3d is 0.5*(tr(C^2)-tr(C)^2;
		// ie C[0][1]*C[1][0] + C[1][2]*C[2][1] + C[2][0]*C[0][2] - C[0][0]*C[1][1] - C[1][1]*C[2][2] - C[2][2]*C[0][0];
		
		ret =   mElementArray[1]*mElementArray[3] + mElementArray[2]*mElementArray[6] 
		      + mElementArray[5]*mElementArray[7] - mElementArray[0]*mElementArray[4]
		      - mElementArray[4]*mElementArray[8] - mElementArray[8]*mElementArray[0];
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
	assert(mRows==mColumns);
	assert(mRows==3);

	return Determinant();
}	

