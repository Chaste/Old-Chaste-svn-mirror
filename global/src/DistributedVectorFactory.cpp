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

#include "DistributedVectorFactory.hpp"
#include "Exception.hpp"

void DistributedVectorFactory::CalculateOwnership(Vec vec)
{
#ifndef NDEBUG
    if (!mPetscStatusKnown)
    {
        CheckForPetsc();
    }
#endif
    // calculate my range
    PetscInt petsc_lo, petsc_hi;
    VecGetOwnershipRange(vec, &petsc_lo, &petsc_hi);
    mLo = (unsigned)petsc_lo;
    mHi = (unsigned)petsc_hi;
    // vector size
    PetscInt size;
    VecGetSize(vec, &size);
    mGlobalHi = (unsigned) size;        
}

DistributedVectorFactory::DistributedVectorFactory(Vec vec) : mPetscStatusKnown(false)
{
    CalculateOwnership(vec);
}

DistributedVectorFactory::DistributedVectorFactory(unsigned size, PetscInt local) : mPetscStatusKnown(false)
{
#ifndef NDEBUG
    if (!mPetscStatusKnown)
    {
        CheckForPetsc();
    }
#endif
    Vec vec;
    VecCreate(PETSC_COMM_WORLD, &vec);
    VecSetSizes(vec, local, size);
    VecSetFromOptions(vec);
    CalculateOwnership(vec);
    VecDestroy(vec);    
}


void DistributedVectorFactory::CheckForPetsc()
{
    assert(mPetscStatusKnown==false);
    PetscTruth petsc_is_initialised;
    PetscInitialized(&petsc_is_initialised);

    //Tripping this assertion means that PETSc and MPI weren't intialised
    //A unit test should include the global fixture:
    //#include "PetscSetupAndFinalize.hpp"
    assert(petsc_is_initialised);
    mPetscStatusKnown=true;
}


Vec DistributedVectorFactory::CreateVec()
{
    Vec vec;
    VecCreate(PETSC_COMM_WORLD, &vec);
    VecSetSizes(vec, mHi-mLo, mGlobalHi);
    VecSetFromOptions(vec);
    return vec;
}

DistributedVector DistributedVectorFactory::CreateDistributedVector(Vec vec)
{
    //Currently the distributed vector class contains static variabled which we need to set here
    ///\todo Make the variables and methods in the DistributedVector class non static
    DistributedVector::SetProblemSize(vec);
    DistributedVector dist_vector(vec); 
    return dist_vector;
}
