#include "MatrixDouble.hpp"
#include <cassert>
#define TINY 0.00000001
#include <cmath>

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
	delete mElementArray;
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





double& MatrixDouble::operator()(int Row, int Column)
{
	assert(Row > -1);
	assert(Row < mRows);
	assert(Column > -1);
	assert(Column < mColumns);
	int Index = Column + mColumns*Row;
	return mElementArray[Index];
}


int MatrixDouble::Rows( void )
{
	return mRows;
}
int MatrixDouble::Columns( void )
{
	return mColumns;
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

double MatrixDouble::Determinant()
{
	assert( mRows == mColumns );
	assert( mRows > 0 && mRows < 4);
	switch( mRows )
	{
		case 1 :
			return mElementArray[0];
			break;
		case 2 :
			return mElementArray[0]*mElementArray[3]-mElementArray[1]*mElementArray[2];
			break;
		case 3 :
			return mElementArray[0]*(mElementArray[8]*mElementArray[4] - mElementArray[5]*mElementArray[7])
             - mElementArray[3]*(mElementArray[1]*mElementArray[8] - mElementArray[2]*mElementArray[7])
             + mElementArray[6]*(mElementArray[1]*mElementArray[5] - mElementArray[2]*mElementArray[4]);
	}
}

MatrixDouble MatrixDouble::Inverse( void )
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
