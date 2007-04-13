#ifndef VECTORSTRIPE_HPP_
#define VECTORSTRIPE_HPP_

#include <vector>
#include <petscvec.h>
#include <iostream>

#include "ParallelVector.hpp"
#include "ParallelIterator.hpp"

class VectorStripe
{
private:
    double *mpVec;
    unsigned mStripe;
    unsigned mStride;
public:
    double& operator[](ParallelIterator index);
    VectorStripe(ParallelVector vector, unsigned stripe);
    
};

#endif /*VECTORSTRIPE_HPP_*/
