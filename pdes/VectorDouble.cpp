#include "VectorDouble.hpp"
#include <cassert>

VectorDouble::VectorDouble(int Size)
{
	assert(Size > 0);
	mSize = Size;
	mElementArray = new double [mSize];
	
	for (int i=0; i<mSize; i++)
	{
		mElementArray[i] = 0.0;
	}
}




VectorDouble::VectorDouble(const VectorDouble& rOtherVector)
{
	mSize = rOtherVector.mSize;
	mElementArray = new double [mSize];
	
	for (int i=0; i<mSize; i++)
	{
		mElementArray[i] = rOtherVector.mElementArray[i];
	}
}







VectorDouble::~VectorDouble()
{
	delete mElementArray;
}







VectorDouble& VectorDouble::operator=(const VectorDouble& rOtherVector)
{
	assert(mSize == rOtherVector.mSize);
	for (int i=0; i<mSize; i++)
	{
		mElementArray[i] = rOtherVector.mElementArray[i];
	}
	return *this;
}





double& VectorDouble::operator()(int Entry)
{
	assert(Entry > -1);
	assert(Entry < mSize);
	return mElementArray[Entry];
}

int VectorDouble::Size(void)
{
	return mSize;
}
