/** Linear System implementation.
*
*
*/
#include "LinearSystem.hpp"
#include "AbstractLinearSolver.hpp"
#include "petscvec.h"
#include "petscksp.h"
#include <iostream>

LinearSystem::LinearSystem(int lhsVectorSize)
{
	//VecCreate(PETSC_COMM_WORLD, &mLhsVector);
    //VecSetSizes(mLhsVector, PETSC_DECIDE, lhsVectorSize);
    //VecSetType(mLhsVector, VECSEQ);
    //VecSetType(mLhsVector, VECMPI);
    //VecSetFromOptions(mLhsVector);
    
    VecCreate(PETSC_COMM_WORLD, &mRhsVector);
    VecSetSizes(mRhsVector, PETSC_DECIDE, lhsVectorSize);
    //VecSetType(mRhsVector, VECSEQ);
    //VecSetType(mRhsVector, VECMPI);
	VecSetFromOptions(mRhsVector);
    
    MatCreate(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,lhsVectorSize,lhsVectorSize,&mLhsMatrix);
    //MatSetType(mLhsMatrix, MATSEQDENSE);
    MatSetType(mLhsMatrix, MATMPIDENSE);
    MatSetFromOptions(mLhsMatrix);
    
    ///\todo: Sparsify matrices (and get the allocation rule correct). 
    //MatSetType(mLhsMatrix, MATMPIAIJ);
    //MatMPIAIJSetPreallocation(mLhsMatrix, 5, PETSC_NULL, 5, PETSC_NULL);
    
	mSize = lhsVectorSize;
    
    VecGetOwnershipRange(mRhsVector, &mOwnershipRangeLo, &mOwnershipRangeHi);
}

LinearSystem::~LinearSystem()
{
	VecDestroy(mRhsVector);
	MatDestroy(mLhsMatrix);
}

//bool LinearSystem::IsMatrixEqualTo(Mat testMatrix)
//{
//    PetscTruth testValue;
//    MatEqual(mLhsMatrix,testMatrix,&testValue);
//    
//    return(testValue == PETSC_TRUE);       
//}
//
//bool LinearSystem::IsRhsVectorEqualTo(Vec testVector)
//{
//   PetscTruth testValue;
//   VecEqual(mRhsVector,testVector, &testValue);
//   
//   return(testValue == PETSC_TRUE);   
//}
void LinearSystem::SetMatrixElement(int row, int col, double value)
{
    if(row >= mOwnershipRangeLo && row < mOwnershipRangeHi)
    {
 		MatSetValue(mLhsMatrix, row, col, value, INSERT_VALUES);
    }
}

void LinearSystem::AddToMatrixElement(int row, int col, double value)
{
    if(row >= mOwnershipRangeLo && row < mOwnershipRangeHi)
    {
 	   MatSetValue(mLhsMatrix, row, col, value, ADD_VALUES);
    }
}

void LinearSystem::AssembleFinalMatrix()
{
    MatAssemblyBegin(mLhsMatrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mLhsMatrix, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(mRhsVector);
    VecAssemblyEnd(mRhsVector);
}

void LinearSystem::AssembleIntermediateMatrix()
{
    MatAssemblyBegin(mLhsMatrix, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(mLhsMatrix, MAT_FLUSH_ASSEMBLY);
    VecAssemblyBegin(mRhsVector);
    VecAssemblyEnd(mRhsVector);
}



void LinearSystem::SetRhsVectorElement(int row, double value)
{
	if(row >= mOwnershipRangeLo && row < mOwnershipRangeHi)
    {    
 	    VecSetValues(mRhsVector, 1, &row, &value, INSERT_VALUES);
    }
}

void LinearSystem::AddToRhsVectorElement(int row, double value)
{
    if(row >= mOwnershipRangeLo && row < mOwnershipRangeHi)    
 	{    
 	   VecSetValues(mRhsVector, 1, &row, &value, ADD_VALUES);
    }
}

void LinearSystem::DisplayMatrix()
{
     MatView(mLhsMatrix,PETSC_VIEWER_STDOUT_WORLD);
}

void LinearSystem::DisplayRhs()
{
     VecView(mRhsVector,PETSC_VIEWER_STDOUT_WORLD);
}

void LinearSystem::SetMatrixRow(int row, double value)
{
    if(row >= mOwnershipRangeLo && row < mOwnershipRangeHi)
    {    
 		int rows, cols;
    	MatGetSize(mLhsMatrix, &rows, &cols);
    	for (int i=0; i<cols; i++)
    	{
        	this->SetMatrixElement(row, i, value);
    	}
    }
}

Vec LinearSystem::Solve(AbstractLinearSolver *solver)
{
    return solver->Solve(mLhsMatrix,mRhsVector);
}

int LinearSystem::GetSize()
{
	return mSize;
}


// this function has been grandfathered because it deadset wrong
/*double LinearSystem::GetMatrixElement(int row, int col)
{
	int row_as_array[1]; row_as_array[0] = row;
	int col_as_array[1]; col_as_array[0] = col;

	double ret_array[0];
	
	MatGetValues(mLhsMatrix, 1, row_as_array, 1, col_as_array, ret_array);

	return ret_array[0];
}
*/
