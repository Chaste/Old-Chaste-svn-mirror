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

#include "LinearSystem.hpp"
#include "PetscException.hpp"
#include <iostream>
#include "OutputFileHandler.hpp"
#include "PetscTools.hpp"
#include <cassert>
#include "HeartEventHandler.hpp"
#include "Timer.hpp"


///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

LinearSystem::LinearSystem(PetscInt lhsVectorSize, MatType matType)
   :mSize(lhsVectorSize),
    mMatNullSpace(NULL),
    mDestroyMatAndVec(true),
    mKspIsSetup(false),
    mNonZerosUsed(0.0),
    mMatrixIsConstant(false),
    mTolerance(1e-6),
    mUseAbsoluteTolerance(false),
    mDirichletBoundaryConditionsVector(NULL),
    mpBlockDiagonalPC(NULL)
{
    SetupVectorAndMatrix(matType);

    /// \todo: if we create a linear system object outside a cardiac assembler, these are gonna
    /// be the default solver and preconditioner. Not consitent with ChasteDefaults.xml though...
    mKspType = "gmres";
    mPcType = "jacobi";

#ifdef TRACE_KSP
    mNumSolves = 0;
    mTotalNumIterations = 0;
    mMaxNumIterations = 0;
#endif
}

LinearSystem::LinearSystem(PetscInt lhsVectorSize, Mat lhsMatrix, Vec rhsVector, MatType matType)
   :mSize(lhsVectorSize),
    mMatNullSpace(NULL),
    mDestroyMatAndVec(true),
    mKspIsSetup(false),
    mNonZerosUsed(0.0),
    mMatrixIsConstant(false),
    mTolerance(1e-6),
    mUseAbsoluteTolerance(false),
    mDirichletBoundaryConditionsVector(NULL),
    mpBlockDiagonalPC(NULL)
{
    // Conveniently, PETSc Mats and Vecs are actually pointers
    mLhsMatrix = lhsMatrix;
    mRhsVector = rhsVector;

    VecGetOwnershipRange(mRhsVector, &mOwnershipRangeLo, &mOwnershipRangeHi);

#ifdef TRACE_KSP
    mNumSolves = 0;
    mTotalNumIterations = 0;
    mMaxNumIterations = 0;
#endif
}

LinearSystem::LinearSystem(Vec templateVector)
   :mMatNullSpace(NULL),
    mDestroyMatAndVec(true),
    mKspIsSetup(false),
    mMatrixIsConstant(false),
    mTolerance(1e-6),
    mUseAbsoluteTolerance(false),
    mDirichletBoundaryConditionsVector(NULL),
    mpBlockDiagonalPC(NULL)
{
    VecDuplicate(templateVector, &mRhsVector);
    VecGetSize(mRhsVector, &mSize);
    VecGetOwnershipRange(mRhsVector, &mOwnershipRangeLo, &mOwnershipRangeHi);
    PetscInt local_size = mOwnershipRangeHi - mOwnershipRangeLo;

    PetscTools::SetupMat(mLhsMatrix, mSize, mSize, (MatType) MATMPIAIJ, local_size, local_size);

    /// \todo: if we create a linear system object outside a cardiac assembler, these are gonna
    /// be the default solver and preconditioner. Not consitent with ChasteDefaults.xml though...
    mKspType = "gmres";
    mPcType = "jacobi";

#ifdef TRACE_KSP
    mNumSolves = 0;
    mTotalNumIterations = 0;
    mMaxNumIterations = 0;
#endif
}

LinearSystem::LinearSystem(Vec residualVector, Mat jacobianMatrix)
    :mMatNullSpace(NULL),
    mDestroyMatAndVec(false),
    mKspIsSetup(false),
    mMatrixIsConstant(false),
    mTolerance(1e-6),
    mUseAbsoluteTolerance(false),
    mDirichletBoundaryConditionsVector(NULL),
    mpBlockDiagonalPC(NULL)
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

    /// \todo: if we create a linear system object outside a cardiac assembler, these are gonna
    /// be the default solver and preconditioner. Not consitent with ChasteDefaults.xml though...
    mKspType = "gmres";
    mPcType = "jacobi";

#ifdef TRACE_KSP
    mNumSolves = 0;
    mTotalNumIterations = 0;
    mMaxNumIterations = 0;
#endif
}

LinearSystem::~LinearSystem()
{
    if (mpBlockDiagonalPC)
    {
        delete mpBlockDiagonalPC;
    }

    if (mDestroyMatAndVec)
    {
        VecDestroy(mRhsVector);
        MatDestroy(mLhsMatrix);
    }

    if (mMatNullSpace)
    {
        MatNullSpaceDestroy(mMatNullSpace);
    }

    if (mKspIsSetup)
    {
        KSPDestroy(mKspSolver);
    }

    if (mDirichletBoundaryConditionsVector)
    {
        ///\todo Never tested in linalg component
        VecDestroy(mDirichletBoundaryConditionsVector);
    }

#ifdef TRACE_KSP
    if (mNumSolves > 0)
    {
        double ave_num_iterations = mTotalNumIterations/(double)mNumSolves;

        std::cout << std::endl << "KSP iterations report:" << std::endl;
        std::cout << "mNumSolves" << "\t" << "mTotalNumIterations" << "\t" << "mMaxNumIterations" << "\t" << "mAveNumIterations" << std::endl;
        std::cout << mNumSolves << "\t" << mTotalNumIterations << "\t" << mMaxNumIterations << "\t" << ave_num_iterations << std::endl;
    }
#endif

}

void LinearSystem::SetupVectorAndMatrix(MatType matType)
{
    VecCreate(PETSC_COMM_WORLD, &mRhsVector);
    VecSetSizes(mRhsVector, PETSC_DECIDE, mSize);
    VecSetFromOptions(mRhsVector);

    PetscTools::SetupMat(mLhsMatrix, mSize, mSize, matType);

    VecGetOwnershipRange(mRhsVector, &mOwnershipRangeLo, &mOwnershipRangeHi);
}

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

#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
    IS is;
    ISCreateGeneral(PETSC_COMM_WORLD,1,&row,&is);
    MatZeroRows(mLhsMatrix, is, &diag_zero);
    ISDestroy(is);
#else

    MatZeroRows(mLhsMatrix, 1, &row, diag_zero);
#endif

}

void LinearSystem::ZeroMatrixColumn(PetscInt col)
{
    MatAssemblyBegin(mLhsMatrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mLhsMatrix, MAT_FINAL_ASSEMBLY);

//#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
   // hello Joe.
   // Hello
//#else
    // determine which rows in this column are non-zero (and
    // therefore need to be zeroed)
    std::vector<unsigned> nonzero_rows;
    for (PetscInt row = mOwnershipRangeLo; row < mOwnershipRangeHi; row++)
    {
        if (GetMatrixElement(row, col) != 0.0)
        {
            nonzero_rows.push_back(row);
        }
    }

    // set those rows to be zero by calling MatSetValues
    unsigned size = nonzero_rows.size();
    PetscInt* rows = new PetscInt[size];
    PetscInt cols[1];
    double* zeros = new double[size];

    cols[0] = col;

    for (unsigned i=0; i<size; i++)
    {
        rows[i] = nonzero_rows[i];
        zeros[i] = 0.0;
    }

    MatSetValues(mLhsMatrix, size, rows, 1, cols, zeros, INSERT_VALUES);
    delete [] rows;
    delete [] zeros;
//#endif
}

void LinearSystem::ZeroRhsVector()
{
    double* p_rhs_vector_array;
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

unsigned LinearSystem::GetSize() const
{
    return (unsigned) mSize;
}

void LinearSystem::SetNullBasis(Vec nullBasis[], unsigned numberOfBases)
{
#ifndef NDEBUG
    // Check all the vectors of the base are normal
    for (unsigned vec_index=0; vec_index<numberOfBases; vec_index++)
    {
        PetscReal l2_norm;
        VecNorm(nullBasis[vec_index], NORM_2, &l2_norm);
        if (fabs(l2_norm-1.0) > 1e-08)
        {
            EXCEPTION("One of the vectors in the null space is not normal");
        }
    }

    // Check all the vectors of the base are orthogonal
    for (unsigned vec_index=1; vec_index<numberOfBases; vec_index++)
    {
        // The strategy is to check the (i-1)-th vector against vectors from i to n with VecMDot. This should be the most efficient way of doing it.
        unsigned num_vectors_ahead = numberOfBases-vec_index;
        PetscScalar dot_products[num_vectors_ahead];
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
        VecMDot(num_vectors_ahead, nullBasis[vec_index-1], &nullBasis[vec_index], dot_products);
#else
        VecMDot(nullBasis[vec_index-1], num_vectors_ahead, &nullBasis[vec_index], dot_products);
#endif
        for (unsigned index=0; index<num_vectors_ahead; index++)
        {
            if (fabs(dot_products[index]) > 1e-08 )
            {
                EXCEPTION("The null space is not orthogonal.");
            }
        }

    }

#endif

    PETSCEXCEPT( MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_FALSE, numberOfBases, nullBasis, &mMatNullSpace) );
}

void LinearSystem::GetOwnershipRange(PetscInt& lo, PetscInt& hi)
{
    lo = mOwnershipRangeLo;
    hi = mOwnershipRangeHi;
}

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

double LinearSystem::GetRhsVectorElement(PetscInt row)
{
    assert(mOwnershipRangeLo <= row && row < mOwnershipRangeHi);

    double* p_rhs_vector;
    PetscInt local_index=row-mOwnershipRangeLo;
    VecGetArray(mRhsVector, &p_rhs_vector);
    double answer=p_rhs_vector[local_index];
    VecRestoreArray(mRhsVector, &p_rhs_vector);

    return answer;
}

unsigned LinearSystem::GetNumIterations() const
{
    assert(this->mKspIsSetup);

    PetscInt num_its;
    KSPGetIterationNumber(this->mKspSolver, &num_its);

    return (unsigned) num_its;
}


Vec& LinearSystem::rGetRhsVector()
{
    return mRhsVector;
}

Vec LinearSystem::GetRhsVector() const
{
    return mRhsVector;
}

Mat& LinearSystem::rGetLhsMatrix()
{
    return mLhsMatrix;
}

Mat LinearSystem::GetLhsMatrix() const
{
    return mLhsMatrix;
}

Vec& LinearSystem::rGetDirichletBoundaryConditionsVector()
{
    return mDirichletBoundaryConditionsVector;
}

void LinearSystem::SetMatrixIsSymmetric(bool isSymmetric)
{
    if (isSymmetric)
    {
#if (PETSC_VERSION_MAJOR == 3) //PETSc 3.x.x
        MatSetOption(mLhsMatrix, MAT_SYMMETRIC, PETSC_TRUE);
        MatSetOption(mLhsMatrix, MAT_SYMMETRY_ETERNAL, PETSC_TRUE);
#else
        MatSetOption(mLhsMatrix, MAT_SYMMETRIC);
        MatSetOption(mLhsMatrix, MAT_SYMMETRY_ETERNAL);
#endif
    }
    else
    {
#if (PETSC_VERSION_MAJOR == 3) //PETSc 3.x.x
        MatSetOption(mLhsMatrix, MAT_SYMMETRIC, PETSC_FALSE);
        MatSetOption(mLhsMatrix, MAT_STRUCTURALLY_SYMMETRIC, PETSC_FALSE);
        MatSetOption(mLhsMatrix, MAT_SYMMETRY_ETERNAL, PETSC_FALSE);
#else
        MatSetOption(mLhsMatrix, MAT_NOT_SYMMETRIC);
        MatSetOption(mLhsMatrix, MAT_NOT_STRUCTURALLY_SYMMETRIC);
        MatSetOption(mLhsMatrix, MAT_NOT_SYMMETRY_ETERNAL);
#endif
    }
}

void LinearSystem::SetMatrixIsConstant(bool matrixIsConstant)
{
    mMatrixIsConstant=matrixIsConstant;
}

void LinearSystem::SetRelativeTolerance(double relativeTolerance)
{
    mTolerance=relativeTolerance;
    mUseAbsoluteTolerance=false;
    if (mKspIsSetup)
    {
        KSPSetTolerances(mKspSolver, mTolerance, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    }
}

void LinearSystem::SetAbsoluteTolerance(double absoluteTolerance)
{
    mTolerance=absoluteTolerance;
    mUseAbsoluteTolerance=true;
    if (mKspIsSetup)
    {
        KSPSetTolerances(mKspSolver, DBL_EPSILON, mTolerance, PETSC_DEFAULT, PETSC_DEFAULT);
    }
}

void LinearSystem::SetKspType(const char *kspType)
{
    mKspType = kspType;
    if (mKspIsSetup)
    {
        KSPSetType(mKspSolver, kspType);
        KSPSetFromOptions(mKspSolver);
    }
}

void LinearSystem::SetPcType(const char* pcType)
{
    mPcType=pcType;
    if (mKspIsSetup)
    {
        if (mPcType == "blockdiagonal")
        {
            if (mpBlockDiagonalPC)
            {
                // If the preconditioner has been set to "blockdiagonal" before, we need to free the pointer.
                delete mpBlockDiagonalPC;
            }
            mpBlockDiagonalPC = new PCBlockDiagonal(mKspSolver);
        }
        else if (mPcType == "ldufactorisation")
        {
            if (mpBlockDiagonalPC)
            {
                // If the preconditioner has been set to "blockdiagonal" before, we need to free the pointer.
                delete mpLDUFactorisationPC;
            }
            mpLDUFactorisationPC = new PCLDUFactorisation(mKspSolver);
        }
        else
        {
            PC prec;
            KSPGetPC(mKspSolver, &prec);
            PCSetType(prec, pcType);
        }
        KSPSetFromOptions(mKspSolver);
    }
}

Vec LinearSystem::Solve(Vec lhsGuess)
{
    /* The following lines are very useful for debugging
     *    MatView(mLhsMatrix,    PETSC_VIEWER_STDOUT_WORLD);
     *    VecView(mRhsVector,    PETSC_VIEWER_STDOUT_WORLD);
     */
    //Double check that the non-zero pattern hasn't changed
    MatInfo mat_info;
    MatGetInfo(mLhsMatrix, MAT_GLOBAL_SUM, &mat_info);

    if (!mKspIsSetup)
    {
        HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
        mNonZerosUsed=mat_info.nz_used;
        //MatNorm(mLhsMatrix, NORM_FROBENIUS, &mMatrixNorm);
        PC prec; //Type of pre-conditioner

        KSPCreate(PETSC_COMM_WORLD, &mKspSolver);
        //See
        //http://www-unix.mcs.anl.gov/petsc/petsc-2/snapshots/petsc-current/docs/manualpages/KSP/KSPSetOperators.html
        //The preconditioner flag (last argument) in the following calls says
        //how to reuse the preconditioner on subsequent iterations
        if (mMatrixIsConstant)
        {
            KSPSetOperators(mKspSolver, mLhsMatrix, mLhsMatrix, SAME_PRECONDITIONER);
        }
        else
        {
            KSPSetOperators(mKspSolver, mLhsMatrix, mLhsMatrix, SAME_NONZERO_PATTERN);
        }

        // Set either absolute or relative tolerance of the KSP solver.
        // The default is to use relative tolerance (1e-6)
        if (mUseAbsoluteTolerance)
        {
            KSPSetTolerances(mKspSolver, DBL_EPSILON, mTolerance, PETSC_DEFAULT, PETSC_DEFAULT);
        }
        else
        {
            KSPSetTolerances(mKspSolver, mTolerance, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
        }

        // set ksp and pc types
        KSPSetType(mKspSolver, mKspType.c_str());
        KSPGetPC(mKspSolver, &prec);


        // Turn off pre-conditioning if the system size is very small
        if (mSize <= 4)
        {
            PCSetType(prec, PCNONE);
        }
        else
        {
            if (mPcType == "blockdiagonal")
            {
                mpBlockDiagonalPC = new PCBlockDiagonal(mKspSolver);
            }
            else if (mPcType == "ldufactorisation")
            {
                if (mpBlockDiagonalPC)
                {
                    // If the preconditioner has been set to "blockdiagonal" before, we need to free the pointer.
                    delete mpLDUFactorisationPC;
                }
                mpLDUFactorisationPC = new PCLDUFactorisation(mKspSolver);
            }
            else
            {
                PCSetType(prec, mPcType.c_str());
            }
        }

        if (mMatNullSpace)
        {
            ///\todo never tested in linalg component
            PETSCEXCEPT( KSPSetNullSpace(mKspSolver, mMatNullSpace) );
        }

        if (lhsGuess)
        {
            //Assume that the user of this method will always be kind enough
            //to give us a reasonable guess.
            KSPSetInitialGuessNonzero(mKspSolver,PETSC_TRUE);
        }

        KSPSetFromOptions(mKspSolver);
        KSPSetUp(mKspSolver);

        mKspIsSetup = true;

        HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);
    }
    else
    {
        #define COVERAGE_IGNORE
        if (mNonZerosUsed!=mat_info.nz_used)
        {
            EXCEPTION("LinearSystem doesn't allow the non-zero pattern of a matrix to change. (I think you changed it).");
        }
//        PetscScalar norm;
//        MatNorm(mLhsMatrix, NORM_FROBENIUS, &norm);
//        if (fabs(norm - mMatrixNorm) > 0)
//        {
//            EXCEPTION("LinearSystem doesn't allow the matrix norm to change");
//        }
        #undef COVERAGE_IGNORE
    }

    // Create solution vector
    ///\todo Should it be compulsory for the caller to supply this and manage the memory?
    HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
    Vec lhs_vector;
    VecDuplicate(mRhsVector, &lhs_vector);//Sets the same size (doesn't copy)
    if (lhsGuess)
    {
        VecCopy(lhsGuess, lhs_vector);
    }
    HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);
//    //Double check that the mRhsVector contains sensible values
//    double* p_rhs,* p_guess;
//    VecGetArray(mRhsVector, &p_rhs);
//    VecGetArray(lhsGuess, &p_guess);
//    for (int global_index=mOwnershipRangeLo; global_index<mOwnershipRangeHi; global_index++)
//    {
//        int local_index = global_index - mOwnershipRangeLo;
//        assert(!std::isnan(p_rhs[local_index]));
//        assert(!std::isnan(p_guess[local_index]));
//        if (p_rhs[local_index] != p_rhs[local_index])
//        {
//            std::cout << "********* PETSc nan\n";
//            assert(0);
//        }
//    }
//    std::cout << "b[0] = " << p_rhs[0] << ", guess[0] = " << p_guess[0] << "\n";
//    VecRestoreArray(mRhsVector, &p_rhs);
//    VecRestoreArray(lhsGuess, &p_guess);
//    // Test A*guess
//    Vec temp;
//    VecDuplicate(mRhsVector, &temp);
//    MatMult(mLhsMatrix, lhs_vector, temp);
//    double* p_temp;
//    VecGetArray(temp, &p_temp);
//    std::cout << "temp[0] = " << p_temp[0] << "\n";
//    VecRestoreArray(temp, &p_temp);
//    VecDestroy(temp);
//    // Dump the matrix to file
//    PetscViewer viewer;
//    PetscViewerASCIIOpen(PETSC_COMM_WORLD,"mat.output",&viewer);
//    MatView(mLhsMatrix, viewer);
//    PetscViewerFlush(viewer);
//    PetscViewerDestroy(viewer);
//    // Dump the rhs vector to file
//    PetscViewerASCIIOpen(PETSC_COMM_WORLD,"vec.output",&viewer);
//    PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
//    VecView(mRhsVector, viewer);
//    PetscViewerFlush(viewer);
//    PetscViewerDestroy(viewer);

    try
    {
        HeartEventHandler::BeginEvent(HeartEventHandler::SOLVE_LINEAR_SYSTEM);
        //Timer::Reset();
        PETSCEXCEPT(KSPSolve(mKspSolver, mRhsVector, lhs_vector));
        //Timer::Print("Solve time:");
        HeartEventHandler::EndEvent(HeartEventHandler::SOLVE_LINEAR_SYSTEM);

        // Check that solver converged and throw if not
        KSPConvergedReason reason;
        KSPGetConvergedReason(mKspSolver, &reason);
        KSPEXCEPT(reason);

#ifdef TRACE_KSP
        PetscInt num_it;
        KSPGetIterationNumber(mKspSolver, &num_it);
        std::cout << "++ Solve: " << mNumSolves << " NumIterations: " << num_it << std::endl << std::flush;

        mNumSolves++;
        mTotalNumIterations += num_it;
        if ((unsigned) num_it > mMaxNumIterations)
        {
            mMaxNumIterations = num_it;
        }
#endif

    }
    catch (const Exception& e)
    {
        // Destroy solution vector on error to avoid memory leaks
        VecDestroy(lhs_vector);
        throw e;
    }

    return lhs_vector;
}
