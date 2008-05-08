/*

Copyright (C) University of Oxford, 2008

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


#include "DistributedVector.hpp"

#include <vector>
#include <petscvec.h>
#include <iostream>
#include <assert.h>

unsigned DistributedVector::mLo=0;
unsigned DistributedVector::mHi=0;
unsigned DistributedVector::mGlobalHi=0;

void DistributedVector::SetProblemSizePerProcessor(unsigned size, PetscInt local)
{
    Vec vec;
    VecCreate(PETSC_COMM_WORLD, &vec);
    VecSetSizes(vec, local, size);
    VecSetFromOptions(vec);
    SetProblemSize(vec);
    VecDestroy(vec);
}

void DistributedVector::SetProblemSize(unsigned size)
{
    SetProblemSizePerProcessor(size, PETSC_DECIDE);
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

unsigned DistributedVector::GetProblemSize()
{
    return mGlobalHi;
}

bool DistributedVector::IsGlobalIndexLocal(unsigned globalIndex)
{
    return (mLo <= globalIndex && globalIndex < mHi);
}

Vec DistributedVector::CreateVec()
{
    Vec vec;
    VecCreate(PETSC_COMM_WORLD, &vec);
    VecSetSizes(vec, mHi-mLo, mGlobalHi);
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

double& DistributedVector::operator[](unsigned globalIndex) throw (DistributedVectorException)
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

double& DistributedVector::operator[](Iterator index) throw (DistributedVectorException)
{
    assert(mStride==1);
    return mpVec[index.Local];
}

