#include "VectorPortion.hpp"

#include <vector>
#include <petscvec.h>
#include <iostream>

VectorPortion::VectorPortion(Vec vec) : mVec(vec)
{
    int lo, hi;
    VecGetOwnershipRange(vec,&lo,&hi);
    mLo = (unsigned) lo;
    mHi = (unsigned) hi;
    VecGetArray(vec, &mpVec);
}

void VectorPortion::Restore()
{
        VecRestoreArray(mVec, &mpVec);
        VecAssemblyBegin(mVec);
        VecAssemblyEnd(mVec);
}

VectorPortion::Iterator VectorPortion::Begin()
{
    Iterator index(this);
    index.Local=0;
    index.Global=mLo;
    return index;
}

VectorPortion::Iterator VectorPortion::End()
{
    Iterator index(this);
    index.Local=mHi-mLo;
    index.Global=mHi;
    return index;
}
