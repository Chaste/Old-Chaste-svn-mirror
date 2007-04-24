#include "DistributedVector.hpp"

#include <vector>
#include <petscvec.h>
#include <iostream>
#include <assert.h>

unsigned DistributedVector::mLo=0;
unsigned DistributedVector::mHi=0;
unsigned DistributedVector::mGlobalHi=0;
void DistributedVector::SetProblemSize(unsigned size)
{
        Vec vec;
        VecCreate(PETSC_COMM_WORLD, &vec);
        VecSetSizes(vec, PETSC_DECIDE, size);
        VecSetFromOptions(vec);
        SetProblemSize(vec);
        VecDestroy(vec);
}

void DistributedVector::SetProblemSize(Vec vec)
{
        // calculate my range
        PetscInt petsc_lo, petsc_hi;
        VecGetOwnershipRange(vec,&petsc_lo,&petsc_hi);
        mLo=(unsigned)petsc_lo;
        mHi=(unsigned)petsc_hi;
        // vector size
        PetscInt size;
        VecGetSize(vec, &size);
        mGlobalHi = (unsigned) size;
}

bool DistributedVector::IsGlobalIndexLocal(unsigned globalIndex)
{
    return (mLo <= globalIndex && globalIndex < mHi);
}

Vec DistributedVector::CreateVec()
{
    Vec vec;
    VecCreate(PETSC_COMM_WORLD, &vec);
    VecSetSizes(vec, PETSC_DECIDE, mGlobalHi);
    VecSetFromOptions(vec);
    return vec;
}

Vec DistributedVector::CreateVec(unsigned stride)
{
    Vec vec;
    VecCreateMPI(PETSC_COMM_WORLD, stride*(mHi-mLo) , stride*mGlobalHi, &vec);
    return vec;
}

DistributedVector::DistributedVector(Vec vec) : mVec(vec)
{
    VecGetArray(vec, &mpVec);
    PetscInt size;
    VecGetSize(vec, &size);
    mStride = (unsigned) size / mGlobalHi;
    assert ((mStride * mGlobalHi) == (unsigned)size);
}

double& DistributedVector::operator[](unsigned globalIndex)
{
    assert(mStride==1);
    if (mLo<=globalIndex && globalIndex <mHi)
    {
        return mpVec[globalIndex - mLo];  
    }
    throw DistributedVectorException();
}

void DistributedVector::Restore()
{
        VecRestoreArray(mVec, &mpVec);
        VecAssemblyBegin(mVec);
        VecAssemblyEnd(mVec);
}

bool DistributedVector::Iterator::operator==(const Iterator& other)
{
    return(Global == other.Global);
}
        
bool DistributedVector::Iterator::operator!=(const Iterator& other)
{
   return(Global != other.Global);
}        

DistributedVector::Iterator& DistributedVector::Iterator::operator++()
{
    Local++;
    Global++;
    return(*this);
}

DistributedVector::Iterator DistributedVector::Begin()
{
    Iterator index;
    index.Local=0;
    index.Global=mLo;
    return index;
}

DistributedVector::Iterator DistributedVector::End()
{
    Iterator index;
    index.Local=mHi-mLo;
    index.Global=mHi;
    return index;
}

double& DistributedVector::operator[](Iterator index)
{
    assert(mStride==1);
    return mpVec[index.Local];
}

