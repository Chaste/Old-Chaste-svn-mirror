#include "ParallelVector.hpp"

#include <vector>
#include <petscvec.h>
#include <iostream>
#include <assert.h>

unsigned ParallelVector::mLo=0;
unsigned ParallelVector::mHi=0;
unsigned ParallelVector::mGlobalHi=0;
void ParallelVector::SetProblemSize(unsigned size)
{
        Vec vec;
        VecCreate(PETSC_COMM_WORLD, &vec);
        VecSetSizes(vec, PETSC_DECIDE, size);
        VecSetFromOptions(vec);
        SetProblemSize(vec);
}

void ParallelVector::SetProblemSize(Vec vec)
{
        // calculate my range
        PetscInt petsc_lo, petsc_hi;
        VecGetOwnershipRange(vec,&petsc_lo,&petsc_hi);
        mLo=(unsigned)petsc_lo;
        mHi=(unsigned)petsc_hi;
        // vector size
        PetscInt size=VecGetSize(vec);
        mGlobalHi = (unsigned) size;
}

ParallelVector::ParallelVector(Vec vec) : mVec(vec)
{
    VecGetArray(vec, &mpVec);
    PetscInt size=VecGetSize(vec);
    mStride = (unsigned) size / mGlobalHi;
    assert ((mStride * mGlobalHi) == (unsigned)size);
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
