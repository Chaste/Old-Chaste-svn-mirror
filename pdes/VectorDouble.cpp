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


double VectorDouble::dot(const VectorDouble &rOtherVector) const
{
	assert(mSize == rOtherVector.Size());
	double product = 0.0;
	for (int i=0; i<mSize; i++)
	{
		product += mElementArray[i] * rOtherVector(i);
	}
	return product;
}


double& VectorDouble::operator()(int Entry) const
{
	assert(Entry > -1);
	assert(Entry < mSize);
	return mElementArray[Entry];
}

int VectorDouble::Size(void) const
{
	return mSize;
}
