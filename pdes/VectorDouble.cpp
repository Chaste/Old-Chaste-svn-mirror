#include "VectorDouble.hpp"
#include <cassert>
#include <iostream>

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

void VectorDouble::ResetToZero( void )
{
	for (int i=0; i<mSize; i++)
	{
		mElementArray[i]=0.0;
	}
}

VectorDouble VectorDouble::VectorProduct(const VectorDouble& rSomeVector )
{
	assert(mSize==3);
	assert(rSomeVector.Size()==3);
	
	VectorDouble result(3);
	
	double x1=(*this)(0);
	double y1=(*this)(1);
	double z1=(*this)(2);	
	double x2=rSomeVector(0);
	double y2=rSomeVector(1);
	double z2=rSomeVector(2);
	
	result(0)=y1*z2-z1*y2;
	result(1)=z1*x2-x1*z2;
	result(2)=x1*y2-y1*x2;
	
	return result;
}
