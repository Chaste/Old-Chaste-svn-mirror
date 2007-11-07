/** Linear System implementation.
*
*
*/
#include "LinearSystem.hpp"
#include "AbstractLinearSolver.hpp"
#include "PetscException.hpp"
#include <iostream>
#include "OutputFileHandler.hpp"
#include "PetscTools.hpp"
#include <cassert>




LinearSystem::LinearSystem(PetscInt lhsVectorSize)
{
    VecCreate(PETSC_COMM_WORLD, &mRhsVector);
    VecSetSizes(mRhsVector, PETSC_DECIDE, lhsVectorSize);
    VecSetFromOptions(mRhsVector);
    
    PetscTools::SetupMat(mLhsMatrix, lhsVectorSize, lhsVectorSize);
    
    ///\todo: Sparsify matrices - get the allocation rule correct.
    //MatMPIAIJSetPreallocation(mLhsMatrix, 5, PETSC_NULL, 5, PETSC_NULL);
    
    mSize = lhsVectorSize;
    
    VecGetOwnershipRange(mRhsVector, &mOwnershipRangeLo, &mOwnershipRangeHi);
    
    mMatNullSpace = NULL;
    mDestroyPetscObjects = true;
}

/**
 * Create a linear system, where the size is based on the size of a given
 * PETSc vec.
 * The LHS & RHS vectors will be created by duplicating this vector's
 * settings.  This should avoid problems with using VecScatter on
 * bidomain simulation results.
 */
LinearSystem::LinearSystem(Vec templateVector)
{
    VecDuplicate(templateVector, &mRhsVector);
    VecGetSize(mRhsVector, &mSize);
    VecGetOwnershipRange(mRhsVector, &mOwnershipRangeLo, &mOwnershipRangeHi);
    PetscInt local_size = mOwnershipRangeHi - mOwnershipRangeLo;

    PetscTools::SetupMat(mLhsMatrix, mSize, mSize, MATMPIAIJ, local_size, local_size);

    mMatNullSpace = NULL;
    mDestroyPetscObjects = true;
}

/**
 * Create a linear system which wraps the provided PETSc objects so we can
 * access them using our API.  Either of the objects may be NULL, but at
 * least one of them must not be.
 * 
 * Useful for storing residuals and jacobians when solving nonlinear PDEs.
 */
LinearSystem::LinearSystem(Vec residualVector, Mat jacobianMatrix)
{
    assert(residualVector || jacobianMatrix);
    mRhsVector = residualVector;
    mLhsMatrix = jacobianMatrix;
    
    PetscInt mat_size=0, vec_size=0;
    if (mRhsVector)
    {
        VecGetSize(mRhsVector, &vec_size);
        mSize = (unsigned)vec_size;
        VecGetOwnershipRange(mRhsVector, &mOwnershipRangeLo, &mOwnershipRangeHi);
    }
    if (mLhsMatrix)
    {
        PetscInt mat_cols;
        MatGetSize(mLhsMatrix, &mat_size, &mat_cols);
        assert(mat_size == mat_cols);
        mSize = (unsigned)mat_size; 
        MatGetOwnershipRange(mLhsMatrix, &mOwnershipRangeLo, &mOwnershipRangeHi);
    }
    assert(!mRhsVector || !mLhsMatrix || vec_size == mat_size);
    
    mMatNullSpace = NULL;
    mDestroyPetscObjects = false;
}

LinearSystem::~LinearSystem()
{
    if (mDestroyPetscObjects)
    {
        VecDestroy(mRhsVector);
        MatDestroy(mLhsMatrix);
    }
    if (mMatNullSpace)
    {
        MatNullSpaceDestroy(mMatNullSpace);
    }
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
void LinearSystem::SetMatrixElement(PetscInt row, PetscInt col, double value)
{
    if (row >= mOwnershipRangeLo && row < mOwnershipRangeHi)
    {
        MatSetValue(mLhsMatrix, row, col, value, INSERT_VALUES);
    }
}

void LinearSystem::AddToMatrixElement(PetscInt row, PetscInt col, double value)
{
    if (row >= mOwnershipRangeLo && row < mOwnershipRangeHi)
    {
        MatSetValue(mLhsMatrix, row, col, value, ADD_VALUES);
    }
}

void LinearSystem::AddToMatrixElements(PetscInt m, PetscInt idxm[], PetscInt n, PetscInt idxn[], double v[])
{
    MatSetValues(mLhsMatrix, m, idxm, n, idxn, v, ADD_VALUES);
}


void LinearSystem::AssembleFinalLinearSystem()
{
    AssembleFinalLhsMatrix();
    AssembleRhsVector();
}


void LinearSystem::AssembleIntermediateLinearSystem()
{
    AssembleIntermediateLhsMatrix();
    AssembleRhsVector();
}

void LinearSystem::AssembleFinalLhsMatrix()
{
    MatAssemblyBegin(mLhsMatrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mLhsMatrix, MAT_FINAL_ASSEMBLY);
}


void LinearSystem::AssembleIntermediateLhsMatrix()
{
    MatAssemblyBegin(mLhsMatrix, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(mLhsMatrix, MAT_FLUSH_ASSEMBLY);
}

void LinearSystem::AssembleRhsVector()
{
    VecAssemblyBegin(mRhsVector);
    VecAssemblyEnd(mRhsVector);
}



void LinearSystem::SetRhsVectorElement(PetscInt row, double value)
{
    if (row >= mOwnershipRangeLo && row < mOwnershipRangeHi)
    {
        VecSetValues(mRhsVector, 1, &row, &value, INSERT_VALUES);
    }
}

void LinearSystem::AddToRhsVectorElement(PetscInt row, double value)
{
    if (row >= mOwnershipRangeLo && row < mOwnershipRangeHi)
    {
        VecSetValues(mRhsVector, 1, &row, &value, ADD_VALUES);
    }
}

void LinearSystem::AddToRhsVectorElements(PetscInt m_rhs, PetscInt idx_rhs[], double rhs[])
{
    VecSetValues(mRhsVector, m_rhs, idx_rhs, rhs, ADD_VALUES);
}

void LinearSystem::DisplayMatrix()
{
    MatView(mLhsMatrix,PETSC_VIEWER_STDOUT_WORLD);
}

void LinearSystem::DisplayRhs()
{
    VecView(mRhsVector,PETSC_VIEWER_STDOUT_WORLD);
}

void LinearSystem::SetMatrixRow(PetscInt row, double value)
{
    if (row >= mOwnershipRangeLo && row < mOwnershipRangeHi)
    {
        PetscInt rows, cols;
        MatGetSize(mLhsMatrix, &rows, &cols);
        for (PetscInt i=0; i<cols; i++)
        {
            this->SetMatrixElement(row, i, value);
        }
    }
}

void LinearSystem::ZeroMatrixRow(PetscInt row)
{
    MatAssemblyBegin(mLhsMatrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mLhsMatrix, MAT_FINAL_ASSEMBLY);
    double diag_zero=0.0;
    // MatZeroRows allows a non-zero value to be placed on the diagonal
    // diag_zero is the value to put in the diagonal
    
#if (PETSC_VERSION_MINOR == 2) //Old API
    IS is;
    ISCreateGeneral(PETSC_COMM_WORLD,1,&row,&is);
    MatZeroRows(mLhsMatrix, is, &diag_zero);
    ISDestroy(is);
#else
    
    MatZeroRows(mLhsMatrix, 1, &row, diag_zero);
#endif
    
}

void LinearSystem::ZeroRhsVector()
{
    double *p_rhs_vector_array;
    VecGetArray(mRhsVector, &p_rhs_vector_array);
    for (PetscInt local_index=0; local_index<mOwnershipRangeHi - mOwnershipRangeLo; local_index++)
    {
        p_rhs_vector_array[local_index]=0.0;
    }
    VecRestoreArray(mRhsVector, &p_rhs_vector_array);
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

Vec LinearSystem::Solve(AbstractLinearSolver *solver, Vec lhsGuess)
{
    return solver->Solve(mLhsMatrix, mRhsVector, mSize, mMatNullSpace, lhsGuess);
}

unsigned LinearSystem::GetSize()
{
    return (unsigned) mSize;
}

void LinearSystem::SetNullBasis(Vec nullBasis[], unsigned numberOfBases)
{
    PETSCEXCEPT( MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_FALSE, numberOfBases, nullBasis, &mMatNullSpace) );
    
    // uncomment following to test null basis is correct:
    //AssembleIntermediateLinearSystem();
    //AssembleFinalLinearSystem();
    //PETSCEXCEPT( MatNullSpaceTest(mMatNullSpace, mLhsMatrix) );
}

/**
 * Get this process' ownership range of the contents of the system
 */
void LinearSystem::GetOwnershipRange(PetscInt &lo, PetscInt &hi)
{
    lo = mOwnershipRangeLo;
    hi = mOwnershipRangeHi;
}

/**
 * Return an element of the matrix.
 * May only be called for elements you own.
 */
double LinearSystem::GetMatrixElement(PetscInt row, PetscInt col)
{
    assert(mOwnershipRangeLo <= row && row < mOwnershipRangeHi);
    PetscInt row_as_array[1];
    row_as_array[0] = row;
    PetscInt col_as_array[1];
    col_as_array[0] = col;
    
    double ret_array[1];
    
    MatGetValues(mLhsMatrix, 1, row_as_array, 1, col_as_array, ret_array);
    
    return ret_array[0];
}

/**
 * Return an element of the RHS vector.
 * May only be called for elements you own.
 */
double LinearSystem::GetRhsVectorElement(PetscInt row)
{
    assert(mOwnershipRangeLo <= row && row < mOwnershipRangeHi);
    
    double *p_rhs_vector;
    PetscInt local_index=row-mOwnershipRangeLo;
    VecGetArray(mRhsVector, &p_rhs_vector);
    double answer=p_rhs_vector[local_index];
    VecRestoreArray(mRhsVector, &p_rhs_vector);
    
    return answer;
}

/**
 * Get access to the rhs vector directly. Shouldn't generally need to be called.
 */
Vec& LinearSystem::rGetRhsVector()
{
    return mRhsVector;
}

/**
 * Get access to the lhs matrix directly. Shouldn't generally need to be called.
 */
Mat& LinearSystem::rGetLhsMatrix()
{
    return mLhsMatrix;
}



/* BROKEN IN PARALLEL
void LinearSystem::WriteLinearSystem(std::string matFile, std::string rhsVectorFile)
{
    OutputFileHandler output_file_handler("");
    out_stream matrix_file = output_file_handler.OpenOutputFile(matFile);
    out_stream vector_file = output_file_handler.OpenOutputFile(rhsVectorFile);

    for(PetscInt i=0; i<mSize; i++)
    {
        for(PetscInt j=0; j<mSize; j++)
        {
            (*matrix_file) << GetMatrixElement(i,j) << " ";
        }
        (*matrix_file) << "\n";
        (*vector_file) << GetRhsVectorElement(i) << "\n";
    }
};
*/
