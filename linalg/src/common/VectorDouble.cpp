#include "VectorDouble.hpp"
#include "Exception.hpp"
#include <cassert>
#include <iostream>
#include <math.h>

/**
 * Default constructor.
 * 
 * This constructor exists so that we can allocate a std::vector<VectorDouble>
 * with a given number of elements.
 * 
 * Note that we don't know what size this VectorDouble will be, so we can't
 * allocate any memory for it.  Hence we set our size to -1, to ensure that
 * an error will be given if this object is used as though it were a normal
 * VectorDouble.  It shouldn't be, because default elements of the std::vector
 * will be replaced by 'real' VectorDouble objects.
 */
VectorDouble::VectorDouble()
{
	mSize = -1;
}

VectorDouble::VectorDouble(int Size)
{
	assert(Size > 0); 
    
    switch (Size) {
        case 1:
            break;
        case 2:
            break;
        case 3:
            break;
        case 4:
            break;
        default:
            std::cout << "Invalid vector size: " << Size << std::endl;
            throw Exception("Invalid vector size");
    }
    
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
	delete[] mElementArray;
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

double& VectorDouble::operator()(int entry) const
{
	assert(entry > -1);
	assert(entry < mSize);
	return mElementArray[entry];
}

int VectorDouble::Size(void) const
{
	return mSize;
}

void VectorDouble::ResetToZero(void)
{
	for (int i=0; i<mSize; i++)
	{
		mElementArray[i]=0.0;
	}
}

VectorDouble VectorDouble::VectorProduct(const VectorDouble& rSomeVector)
{
	assert(mSize==3);
	assert(rSomeVector.Size()==3);
	
	VectorDouble result(3);
	
	double x1=mElementArray[0];
	double y1=mElementArray[1];
	double z1=mElementArray[2];
	double x2=rSomeVector(0);
	double y2=rSomeVector(1);
	double z2=rSomeVector(2);
	
	result(0)=y1*z2-z1*y2;
	result(1)=z1*x2-x1*z2;
	result(2)=x1*y2-y1*x2;
	
	return result;
}

VectorDouble VectorDouble::operator*(double scalar)
{
	VectorDouble result(mSize);
	for (int i=0; i<mSize; i++)
	{
		result(i)=scalar*mElementArray[i];
	}
	return result;
}

VectorDouble operator*(double scalar, const VectorDouble& rSomeVector)
{
	VectorDouble result(rSomeVector.Size());
	for (int i=0; i<rSomeVector.Size(); i++)
	{
		result(i)=scalar*rSomeVector(i);
	}
	return result;
}

VectorDouble VectorDouble::operator+(const VectorDouble& rSomeVector1)
{
	assert(this->Size()==rSomeVector1.Size());
	VectorDouble result(this->Size());
	for (int i=0; i<this->Size(); i++)
	{
		result(i)=mElementArray[i]+rSomeVector1(i);
	}
	return result;
}

VectorDouble VectorDouble::operator-(const VectorDouble& rSomeVector1)
{
	assert(this->Size()==rSomeVector1.Size());
	VectorDouble result(this->Size());
	for (int i=0; i<this->Size(); i++)
	{
		result(i)=mElementArray[i]-rSomeVector1(i);
	}
	return result;
}

double VectorDouble::L2Norm(void)
{
	return sqrt((*this).dot(*this));
}
