/*

Copyright (C) University of Oxford, 2005-2009

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


#include "ReplicatableVector.hpp"

#include <cassert>

// Private methods

void ReplicatableVector::RemovePetscContext()
{
    if (mToAll != NULL)
    {
        VecScatterDestroy(mToAll);
        mToAll = NULL;
    }

    if (mReplicated != NULL)
    {
        VecDestroy(mReplicated);
        mReplicated = NULL;
    }

    if (mDistributed != NULL)
    {
        VecDestroy(mDistributed);
        mDistributed = NULL;
    }
}

// Constructors & destructors

ReplicatableVector::ReplicatableVector()
    : mToAll(NULL),
      mReplicated(NULL),
      mDistributed(NULL)
{
}

ReplicatableVector::ReplicatableVector(Vec vec)
    : mToAll(NULL),
      mReplicated(NULL),
      mDistributed(NULL)
{
    ReplicatePetscVector(vec);
}

ReplicatableVector::ReplicatableVector(unsigned size)
    : mToAll(NULL),
      mReplicated(NULL),
      mDistributed(NULL)
{
    resize(size);
}

ReplicatableVector::~ReplicatableVector()
{
    RemovePetscContext();
}


// Vector interface methods

unsigned ReplicatableVector::size()
{
    return mData.size();
}

void ReplicatableVector::resize(unsigned size)
{
    // PETSc stuff will be out of date
    RemovePetscContext();
    mData.resize(size);
}

double& ReplicatableVector::operator[](unsigned index)
{
    assert(index < mData.size());
    return mData[index];
}


// The workhorse methods

void ReplicatableVector::Replicate(unsigned lo, unsigned hi)
{
    // Copy information into a PetSC vector
    if (mDistributed==NULL)
    {
        VecCreateMPI(PETSC_COMM_WORLD, hi-lo, this->size(), &mDistributed);
    }

    double* p_distributed;
    VecGetArray(mDistributed, &p_distributed);
    for (unsigned global_index=lo; global_index<hi; global_index++)
    {
        p_distributed[ (global_index-lo) ] = mData[global_index];
    }
    VecAssemblyBegin(mDistributed);
    VecAssemblyEnd(mDistributed);

    // Now do the real replication
    ReplicatePetscVector(mDistributed);
}

void ReplicatableVector::ReplicatePetscVector(Vec vec)
{
    // If the size has changed then we'll need to make a new context
    PetscInt isize;
    VecGetSize(vec, &isize);
    unsigned size=isize;

    if (this->size() != size)
    {
        resize(size);
    }
    if (mReplicated == NULL)
    {
        // This creates mReplicated (the scatter context) and mReplicated (to store values)
        VecScatterCreateToAll(vec, &mToAll, &mReplicated);
    }

    // Replicate the data
//PETSc-3.x.x or PETSc-2.3.3 
#if ( (PETSC_VERSION_MAJOR == 3) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 3)) //2.3.3 or 3.x.x
    VecScatterBegin(mToAll, vec, mReplicated, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd  (mToAll, vec, mReplicated, INSERT_VALUES, SCATTER_FORWARD);
#else
//PETSc-2.3.2 or previous
    VecScatterBegin(vec, mReplicated, INSERT_VALUES, SCATTER_FORWARD, mToAll);
    VecScatterEnd  (vec, mReplicated, INSERT_VALUES, SCATTER_FORWARD, mToAll);
#endif

    // Information is now in mReplicated PETSc vector
    // Copy into mData
    double* p_replicated;
    VecGetArray(mReplicated, &p_replicated);
    for (unsigned i=0; i<size; i++)
    {
        mData[i] = p_replicated[i];
    }
}
