#include "ReplicatableVector.hpp"

#include <vector>
#include <petscvec.h>
#include <iostream>

// Private methods

void ReplicatableVector::RemovePetscContext()
{
    if (mToAll != NULL)
    {
        VecScatterDestroy(mToAll);
        mToAll=NULL;
    }
    
    if (mReplicated != NULL)
    {
        VecDestroy(mReplicated);
        mReplicated=NULL;
    }
    
    if (mDistributed != NULL)
    {
        VecDestroy(mDistributed);
        mDistributed=NULL;
    }
}

// Constructors & destructors

ReplicatableVector::ReplicatableVector()
{
    mToAll=NULL;
    mReplicated=NULL;
    mDistributed=NULL;
}

ReplicatableVector::ReplicatableVector(Vec vec)
{
    mToAll=NULL;
    mReplicated=NULL;
    mDistributed=NULL;
    
    ReplicatePetscVector(vec);
}

ReplicatableVector::ReplicatableVector(unsigned size)
{
    mToAll=NULL;
    mReplicated=NULL;
    mDistributed=NULL;
    resize(size);
}

ReplicatableVector::~ReplicatableVector()
{
    RemovePetscContext();
}


// Vector interface methods

unsigned ReplicatableVector::size(void)
{
    return mData.size();
}

void ReplicatableVector::resize(unsigned size)
{
    //PETSc stuff will be out of date
    RemovePetscContext();
    mData.resize(size);
}

double& ReplicatableVector::operator[](unsigned index)
{
    return mData[index];
}


// The workhorse methods

void ReplicatableVector::Replicate(unsigned lo, unsigned hi)
{
    //Copy information into a PetSC vector
    if (mDistributed==NULL)
    {
        VecCreateMPI(PETSC_COMM_WORLD, hi-lo, this->size(), &mDistributed);
    }
    
    double *p_distributed;
    VecGetArray(mDistributed, &p_distributed);
    for (unsigned global_index=lo; global_index<hi; global_index++)
    {
        p_distributed[ (global_index-lo) ]= mData[global_index];
    }
    VecAssemblyBegin(mDistributed);
    VecAssemblyEnd(mDistributed);
    
    //Now do the real replication
    ReplicatePetscVector(mDistributed);
}

void ReplicatableVector::ReplicatePetscVector(Vec vec)
{
    //If the size has changed then we'll need to make a new context
    PetscInt isize;
    VecGetSize(vec, &isize);
    unsigned size=isize;
    
    if (this->size() != size)
    {
        resize(size);
    }
    if (mReplicated == NULL)
    {
        //This creates mReplicated (the scatter context) and mReplicated (to store values)
        VecScatterCreateToAll(vec, &mToAll, &mReplicated);
    }
    
    //Replicate the data
    VecScatterBegin(vec, mReplicated, INSERT_VALUES, SCATTER_FORWARD, mToAll);
    VecScatterEnd  (vec, mReplicated, INSERT_VALUES, SCATTER_FORWARD, mToAll);
    
    //Information is now in mReplicated PETSc vector
    //Copy into mData
    double *p_replicated;
    VecGetArray(mReplicated, &p_replicated);
    for (unsigned i=0; i<size; i++)
    {
        mData[i]=p_replicated[i];
    }
}
