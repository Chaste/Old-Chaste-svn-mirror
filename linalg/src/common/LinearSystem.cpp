/** Linear System implementation.
*
*
*/
#include "LinearSystem.hpp"
#include "AbstractLinearSolver.hpp"
#include <iostream>

LinearSystem::LinearSystem(int lhsVectorSize)
{
    
    VecCreate(PETSC_COMM_WORLD, &mRhsVector);
    VecSetSizes(mRhsVector, PETSC_DECIDE, lhsVectorSize);
	VecSetFromOptions(mRhsVector);
    
    MatCreate(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,lhsVectorSize,lhsVectorSize,&mLhsMatrix);
    MatSetType(mLhsMatrix, MATMPIAIJ);
    MatSetFromOptions(mLhsMatrix);
    
    ///\todo: Sparsify matrices - get the allocation rule correct. 
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

void LinearSystem::AssembleFinalLinearSystem()
{
    MatAssemblyBegin(mLhsMatrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mLhsMatrix, MAT_FINAL_ASSEMBLY);
    AssembleRhsVector();
}


void LinearSystem::AssembleIntermediateLinearSystem()
{
    MatAssemblyBegin(mLhsMatrix, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(mLhsMatrix, MAT_FLUSH_ASSEMBLY);
    AssembleRhsVector();
}


void LinearSystem::AssembleRhsVector()
{
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


void LinearSystem::ZeroMatrixRow(int row)
{
    MatAssemblyBegin(mLhsMatrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mLhsMatrix, MAT_FINAL_ASSEMBLY);
    IS is;
    ISCreateGeneral(PETSC_COMM_WORLD,1,&row,&is);
    double zero=0.0;
    MatZeroRows(mLhsMatrix, is, &zero);
    ISDestroy(is);
}

void LinearSystem::ZeroRhsVector()
{
    double *p_rhs_vector;
    VecGetArray(mRhsVector, &p_rhs_vector);
    for (int local_index=0; local_index<mOwnershipRangeHi - mOwnershipRangeLo; local_index++)
    {   
        p_rhs_vector[local_index]=0.0;
    }
    VecRestoreArray(mRhsVector, &p_rhs_vector);  
}


void LinearSystem::ZeroLhsMatrix()
{
    MatZeroEntries(mLhsMatrix);
}

void LinearSystem::ZeroLinearSystem()
{
    ZeroRhsVector();
    ZeroLhsMatrix();
}
 
Vec LinearSystem::Solve(AbstractLinearSolver *solver)
{
    return solver->Solve(mLhsMatrix, mRhsVector, mSize);
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
