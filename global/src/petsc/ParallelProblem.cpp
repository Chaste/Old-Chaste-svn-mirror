#include <petscvec.h>

#include "ParallelProblem.hpp"


void ParallelProblem::SetProblemSize(unsigned problemSize)
{
     // create a 10 element petsc vector
    Vec vec;

    VecCreate(PETSC_COMM_WORLD, &vec);
    VecSetSizes(vec, PETSC_DECIDE, problemSize);
    VecSetFromOptions(vec);
    
    // calculate my range
    PetscInt petsc_lo, petsc_hi;
    VecGetOwnershipRange(vec,&petsc_lo,&petsc_hi);
    mLo=(unsigned)petsc_lo;
    mHi=(unsigned)petsc_hi; 
    mProblemSize = problemSize;  
}

unsigned ParallelProblem::GetProblemSize()
{
    return mProblemSize;
}
unsigned ParallelProblem::Size()
{
    return mHi-mLo;
}

ParallelIterator ParallelProblem::Begin()
{
    return 0;
}

ParallelIterator ParallelProblem::End()
{
    return Size();
}

unsigned ParallelProblem::Global(ParallelIterator iterator)
{
    return mLo+iterator;
}

unsigned ParallelProblem::Local(ParallelIterator iterator)
{
    return iterator;
}
