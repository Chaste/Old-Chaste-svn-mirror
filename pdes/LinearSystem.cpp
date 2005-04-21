/** Linear System implementation.
*
*
*
*/
#include "LinearSystem.hpp"
#include "AbstractLinearSolver.hpp"
#include "petscvec.h"
#include <iostream>

LinearSystem::LinearSystem(int lhsVectorSize)
{
    VecCreate(PETSC_COMM_WORLD, &mLhsVector);
    VecSetSizes(mLhsVector, PETSC_DECIDE, lhsVectorSize);
    VecSetType(mLhsVector, VECSEQ);
    
    VecCreate(PETSC_COMM_WORLD, &mRhsVector);
    VecSetSizes(mRhsVector, PETSC_DECIDE, lhsVectorSize);
    VecSetType(mRhsVector, VECSEQ);

    MatCreate(PETSC_COMM_WORLD,lhsVectorSize,lhsVectorSize,PETSC_DETERMINE,PETSC_DETERMINE,&mLhsMatrix);
    MatSetType(mLhsMatrix, MATSEQDENSE);
}

bool LinearSystem::IsMatrixEqualTo(Mat testMatrix)
{
    PetscTruth testValue;
    MatEqual(mLhsMatrix,testMatrix,&testValue);
    
    return(testValue == PETSC_TRUE);       
}

bool LinearSystem::IsRhsVectorEqualTo(Vec testVector)
{
   PetscTruth testValue;
   VecEqual(mRhsVector,testVector, &testValue);
   
   return(testValue == PETSC_TRUE);   
}
void LinearSystem::SetMatrixElement(int row, int col, double value)
{
    MatSetValue(mLhsMatrix, row, col, value, INSERT_VALUES);
}

void LinearSystem::AddToMatrixElement(int row, int col, double value)
{
    MatSetValue(mLhsMatrix, row, col, value, ADD_VALUES);
}

void LinearSystem::AssembleFinalMatrix()
{
    MatAssemblyBegin(mLhsMatrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mLhsMatrix, MAT_FINAL_ASSEMBLY);
}

void LinearSystem::AssembleIntermediateMatrix()
{
    MatAssemblyBegin(mLhsMatrix, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(mLhsMatrix, MAT_FLUSH_ASSEMBLY);
}



void LinearSystem::SetRhsVectorElement(int row, double value)
{
    VecSetValue(mRhsVector, row, value, INSERT_VALUES);
}

void LinearSystem::AddToRhsVectorElement(int row, double value)
{
    VecSetValue(mRhsVector, row, value, ADD_VALUES);
}

void LinearSystem::DisplayMatrix()
{
     MatView(mLhsMatrix,PETSC_VIEWER_STDOUT_SELF);
}

void LinearSystem::SetMatrixRow(int row, double value)
{
    int rows, cols;
    MatGetSize(mLhsMatrix, &rows, &cols);
    for (int i=0; i<cols; i++)
    {
        this->SetMatrixElement(row, i, value);
    }
}

Vec LinearSystem::Solve(AbstractLinearSolver *solver)
{
    return solver->Solve(mLhsMatrix,mRhsVector);
}
