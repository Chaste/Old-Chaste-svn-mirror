#ifndef VECTORSTRIPE_CPP_
#define VECTORSTRIPE_CPP_

#include "VectorStripe.hpp"

VectorStripe::VectorStripe(ParallelVector vector, unsigned stripe)
{
    mpVec=vector.mpVec;
    mStripe = stripe;
    mStride=vector.mStride;
}


double& VectorStripe::operator[](ParallelIterator index)
{
    return mpVec[mStride*index+mStripe];
}


#endif /*VECTORSTRIPE_CPP_*/
