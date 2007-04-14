#include "ParallelVector.hpp"

#include <vector>
#include <petscvec.h>
#include <iostream>

unsigned ParallelVector::mLo=0;
unsigned ParallelVector::mHi=0;

ParallelVector::ParallelVector(Vec vec) : mVec(vec)
{
    int lo, hi;
    VecGetOwnershipRange(vec,&lo,&hi);
    mLo = (unsigned) lo;
    mHi = (unsigned) hi;
    VecGetArray(vec, &mpVec);
}

void ParallelVector::Restore()
{
        VecRestoreArray(mVec, &mpVec);
        VecAssemblyBegin(mVec);
        VecAssemblyEnd(mVec);
}

ParallelVector::Iterator ParallelVector::Begin()
{
    Iterator index;
    index.Local=0;
    index.Global=mLo;
    return index;
}

ParallelVector::Iterator ParallelVector::End()
{
    Iterator index;
    index.Local=mHi-mLo;
    index.Global=mHi;
    return index;
}

double& ParallelVector::operator[](Iterator index)
{
    return mpVec[index.Local];
}
