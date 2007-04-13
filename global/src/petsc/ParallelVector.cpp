#include "ParallelVector.hpp"
#include "GlobalParallelProblem.hpp"

#include <vector>
#include <petscvec.h>
#include <iostream>

ParallelVector::ParallelVector(Vec vec) : mVec(vec)
{
    int lo, hi;
    int size;
    VecGetSize(vec, &size);
    
    
    VecGetOwnershipRange(vec,&lo,&hi);
    VecGetArray(vec, &mpVec);
    mStride = (unsigned)size / (gProblem.End()-gProblem.Begin());
    std::cout << "mStride" << mStride << std::endl;
}

void ParallelVector::Restore()
{
        VecRestoreArray(mVec, &mpVec);
        VecAssemblyBegin(mVec);
        VecAssemblyEnd(mVec);
}

double& ParallelVector::operator[](ParallelIterator index)
{
    return mpVec[index];
}
