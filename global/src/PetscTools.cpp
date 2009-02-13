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

#include "PetscTools.hpp"

bool PetscTools::mPetscIsInitialised = false;
unsigned PetscTools::mNumProcessors = 0;
unsigned PetscTools::mRank = 0;

void PetscTools::ResetCache()
{
#ifdef SPECIAL_SERIAL
    mPetscIsInitialised = false;
    mNumProcessors = 1;
    mRank = 0;
#else
    PetscTruth is_there;
    PetscInitialized(&is_there);
    if (is_there)
    {
        mPetscIsInitialised = true;
        
        PetscInt num_procs;
        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);
        mNumProcessors = (unsigned) num_procs;
        
        PetscInt my_rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);
        mRank = (unsigned) my_rank;
    }
    else
    {
        // No PETSc
        mPetscIsInitialised = false;
        mNumProcessors = 1;
        mRank = 0;
    }
#endif
}

//
// Information methods
//

bool PetscTools::IsSequential()
{
    if (mNumProcessors == 0)
    {
        ResetCache();
    }
    return (mNumProcessors == 1);
}

unsigned PetscTools::NumProcs()
{
    if (mNumProcessors == 0)
    {
        ResetCache();
    }
    return mNumProcessors;
}

unsigned PetscTools::GetMyRank()
{
    if (mNumProcessors == 0)
    {
        ResetCache();
    }
    return mRank;
}

bool PetscTools::AmMaster()
{
    if (mNumProcessors == 0)
    {
        ResetCache();
    }
    return (mRank == MASTER_RANK);
}

//
// Little utility methods
//

void PetscTools::Barrier()
{
    if (mNumProcessors == 0)
    {
        ResetCache();
    }
    if (mPetscIsInitialised)
    {
        PetscBarrier(PETSC_NULL);
    }
}

#ifndef SPECIAL_SERIAL

void PetscTools::ReplicateException(bool flag)
    {
        unsigned my_error = (unsigned) flag;
        unsigned anyones_error;
        MPI_Allreduce(&my_error, &anyones_error, 1, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);
        if (flag)
        {
            // Return control to exception thrower
            return;
        }
        if (anyones_error)
        {
            EXCEPTION("Another process threw an exception; bailing out.");
        }
    }

//
// Vector & Matrix creation routines
//

Vec PetscTools::CreateVec(int size)
{
    assert(size>0);
    Vec ret;
    VecCreate(PETSC_COMM_WORLD, &ret);
    VecSetSizes(ret, PETSC_DECIDE, size);
    VecSetFromOptions(ret);
    return ret;
}

Vec PetscTools::CreateVec(int size, double value)
{
    assert(size>0);
    Vec ret = CreateVec(size);

    #if (PETSC_VERSION_MINOR == 2) //Old API
    VecSet(&value, ret);
    #else
    VecSet(ret, value);
    #endif

    VecAssemblyBegin(ret);
    VecAssemblyEnd(ret);
    return ret;
}

Vec PetscTools::CreateVec(std::vector<double> data)
{
    assert(data.size()>0);
    Vec ret = CreateVec(data.size());

    double* p_ret;
    VecGetArray(ret, &p_ret);
    int lo, hi;
    VecGetOwnershipRange(ret, &lo, &hi);

    for (int global_index=lo; global_index < hi; global_index++)
    {
        int local_index = global_index - lo;
        p_ret[local_index] = data[global_index];
    }
    VecRestoreArray(ret, &p_ret);
    VecAssemblyBegin(ret);
    VecAssemblyEnd(ret);

    return ret;
}

void PetscTools::SetupMat(Mat& rMat, int numRows, int numColumns,
                          MatType matType,
                          int numLocalRows,
                          int numLocalColumns,
                          int maxColsPerRow)
{
    assert(numRows>0);
    assert(numColumns>0);

    #if (PETSC_VERSION_MINOR == 2) //Old API
    MatCreate(PETSC_COMM_WORLD,numLocalRows,numLocalColumns,numRows,numColumns,&rMat);
    #else //New API
    MatCreate(PETSC_COMM_WORLD,&rMat);
    MatSetSizes(rMat,numLocalRows,numLocalColumns,numRows,numColumns);
    #endif

    MatSetType(rMat, matType);

    if (strcmp(matType,MATMPIAIJ)==0)
    {
        MatMPIAIJSetPreallocation(rMat, maxColsPerRow, PETSC_NULL, (PetscInt) (maxColsPerRow*0.5), PETSC_NULL);
    }

    MatSetFromOptions(rMat);
}

#endif //SPECIAL_SERIAL
