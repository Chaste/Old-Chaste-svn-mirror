#include "MatrixDouble.hpp"
#include <cassert>

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
