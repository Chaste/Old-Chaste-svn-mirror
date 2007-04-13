#include "ParallelVector.hpp"
#include "GlobalParallelProblem.hpp"

#include <vector>
#include <petscvec.h>
#include <iostream>
#include <assert.h>

ParallelVector::ParallelVector(Vec vec) : mVec(vec)
{
    int lo, hi;
    int vector_size;
    VecGetSize(vec, &vector_size);
    
    
    VecGetOwnershipRange(vec,&lo,&hi);
    VecGetArray(vec, &mpVec);
    mStride = (unsigned)vector_size / (gProblem.GetProblemSize());
    assert(gProblem.GetProblemSize()*mStride == (unsigned)vector_size);
}

void ParallelVector::Restore()
{
        VecRestoreArray(mVec, &mpVec);
        VecAssemblyBegin(mVec);
        VecAssemblyEnd(mVec);
}

double& ParallelVector::operator[](ParallelIterator index)
{
    assert(mStride == 1);
    return mpVec[index];
}
