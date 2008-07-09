/*

Copyright (C) University of Oxford, 2008

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


#ifndef _LINEARSYSTEM_HPP_
#define _LINEARSYSTEM_HPP_

#include "UblasCustomFunctions.hpp"
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>

#include <string>

/**
 * Linear System class. Stores and solves a linear equation of the form Ax=b,
 * where A is a square matrix and x and b are column vectors.
 * The class uses PETSc.
 */
class LinearSystem
{
	friend class TestLinearSystem;

private:
    Mat mLhsMatrix;
    Vec mRhsVector;
    PetscInt mSize;
    /** \todo
     * Verify claim that ownership range for Vec and Mat is same.
     * This should only matter for efficiency if the claim is false.
     */
    PetscInt mOwnershipRangeLo; /*< For parallel code.  Stores lowest index of vectors and lowest row of matrix*/
    PetscInt mOwnershipRangeHi; /*< Stores <b>one more than</b> the highest index stored locally*/

    MatNullSpace mMatNullSpace;

    /** Whether we need to destroy the Petsc matrix and vector in our destructor */
    bool mDestroyMatAndVec;

    KSP mKspSolver;
    bool mKspIsSetup; //Used by Solve method to track whether KSP has been used
    double mNonZerosUsed; //Yes, it really is stored as a double.
    bool mMatrixIsConstant;
    double mTolerance;
    bool mUseAbsoluteTolerance;
    char mKspType[30];
    char mPcType[30];
    
public:
    LinearSystem(PetscInt lhsVectorSize);
    LinearSystem(Vec templateVector);
    LinearSystem(Vec residualVector, Mat jacobianMatrix);
    ~LinearSystem();

//    bool IsMatrixEqualTo(Mat testMatrix);
//    bool IsRhsVectorEqualTo(Vec testVector);
    void SetMatrixElement(PetscInt row, PetscInt col, double value);
    void AddToMatrixElement(PetscInt row, PetscInt col, double value);

    void AssembleFinalLinearSystem();         // Call before solve
    void AssembleIntermediateLinearSystem();  // Should be called before AddToMatrixElement
    void AssembleFinalLhsMatrix();
    void AssembleIntermediateLhsMatrix();
    void AssembleRhsVector();

    void SetMatrixIsSymmetric();
    void SetMatrixIsConstant(bool matrixIsConstant);
    void SetRelativeTolerance(double relativeTolerance);
    void SetAbsoluteTolerance(double absoluteTolerance);
    void SetKspType(const char*);
    void SetPcType(const char*);
    void DisplayMatrix();
    void DisplayRhs();
    void SetMatrixRow(PetscInt row, double value);
    void ZeroMatrixRow(PetscInt row);
    void ZeroLhsMatrix();
    void ZeroRhsVector();
    void ZeroLinearSystem();
    Vec Solve(Vec lhsGuess=NULL);
    void SetRhsVectorElement(PetscInt row, double value);
    void AddToRhsVectorElement(PetscInt row, double value);
    unsigned GetSize();
    void SetNullBasis(Vec nullbasis[], unsigned numberOfBases);
    Vec& rGetRhsVector();
    Mat& rGetLhsMatrix();


    // DEBUGGING CODE:
    void GetOwnershipRange(PetscInt &lo, PetscInt &hi);
    double GetMatrixElement(PetscInt row, PetscInt col);
    double GetRhsVectorElement(PetscInt row);


    /***
     * Add multiple values to the matrix of linear system
     * @param matrixRowAndColIndices mapping from index of the ublas matrix (see param below)
     *  to index of the Petsc matrix of this linear system
     * @param smallMatrix Ublas matrix containing the values to be added
     *
     * N.B. Values which are not local (ie the row is not owned) will be skipped.
     */
    template<size_t MATRIX_SIZE>
    void AddLhsMultipleValues(unsigned* matrixRowAndColIndices, c_matrix<double, MATRIX_SIZE, MATRIX_SIZE>& smallMatrix)
    {
        PetscInt matrix_row_indices[MATRIX_SIZE];
        PetscInt num_rows_owned=0;
        unsigned num_values_owned=0;

        double values[MATRIX_SIZE*MATRIX_SIZE];
        for (unsigned row = 0; row<MATRIX_SIZE; row++)
        {
            PetscInt global_row = matrixRowAndColIndices[row];
            if (global_row >=mOwnershipRangeLo && global_row <mOwnershipRangeHi)
            {
                matrix_row_indices[num_rows_owned++] = global_row;
                for (unsigned col=0; col<MATRIX_SIZE; col++)
                {
                    values[num_values_owned++] = smallMatrix(row,col);
                }
            }
        }

        MatSetValues(mLhsMatrix,
                     num_rows_owned,
                     matrix_row_indices,
                     MATRIX_SIZE,
                     (PetscInt*) matrixRowAndColIndices,
                     values,
                     ADD_VALUES);
    };


    /***
     * Add multiple values to the RHS vector
     * @param vectorIndices mapping from index of the ublas vector (see param below)
     *  to index of the vector of this linear system
     * @param smallVector Ublas vector containing the values to be added
     *
     * N.B. Values which are not local (ie the row is not owned) will be skipped.
     */
    template<size_t VECTOR_SIZE>
    void AddRhsMultipleValues(unsigned* VectorIndices, c_vector<double, VECTOR_SIZE>& smallVector)
    {
        PetscInt indices_owned[VECTOR_SIZE];
        PetscInt num_indices_owned=0;

        double values[VECTOR_SIZE];
        for (unsigned row = 0; row<VECTOR_SIZE; row++)
        {
            PetscInt global_row = VectorIndices[row];
            if (global_row >=mOwnershipRangeLo && global_row <mOwnershipRangeHi)
            {
                indices_owned[num_indices_owned] = global_row;
                values[num_indices_owned] = smallVector(row);
                num_indices_owned++;
            }
        }

        VecSetValues(mRhsVector,
                     num_indices_owned,
                     indices_owned,
                     values,
                     ADD_VALUES);
    }
};

#endif //_LINEARSYSTEM_HPP_
