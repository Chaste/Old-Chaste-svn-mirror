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


#ifndef _LINEARSYSTEM_HPP_
#define _LINEARSYSTEM_HPP_

#include <boost/serialization/access.hpp>
#include "UblasCustomFunctions.hpp"
#include "PetscTools.hpp"
#include "OutputFileHandler.hpp"
#include "PCBlockDiagonal.hpp"
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscviewer.h>

#include <string>
#include <cassert>
 
// Needs to be included last
#include <boost/serialization/export.hpp>

/**
 * Linear System class. Stores and solves a linear equation of the form Ax=b,
 * where A is a square matrix and x and b are column vectors.
 * The class uses PETSc.
 */
class LinearSystem
{
    friend class TestLinearSystem;

private:

    Mat mLhsMatrix;  /**< The left-hand side matrix. */
    Vec mRhsVector;  /**< The right-hand side vector. */
    PetscInt mSize;  /**< The size of the linear system. */

    /** \todo
     * Verify claim that ownership range for Vec and Mat is same.
     * This should only matter for efficiency if the claim is false.
     */
    PetscInt mOwnershipRangeLo; /**< For parallel code.  Stores lowest index of vectors and lowest row of matrix stored locally. */
    PetscInt mOwnershipRangeHi; /**< Stores <b>one more than</b> the highest index stored locally. */

    MatNullSpace mMatNullSpace; /**< PETSc null matrix. */

    /** Whether we need to destroy the PETSc matrix and vector in our destructor */
    bool mDestroyMatAndVec;

    KSP mKspSolver;   /**< The PETSc linear solver object */
    bool mKspIsSetup; /**< Used by Solve method to track whether KSP has been used. */
    double mNonZerosUsed;  /**< Yes, it really is stored as a double. */
    bool mMatrixIsConstant; /**< Whether the matrix is unchanged each time Solve() is called */
    double mTolerance; /**< absolute or relative tolerance of the KSP solver */
    /**
     * Sets either absolute or relative tolerance of the KSP solver.
     * Default is to false
     */
    bool mUseAbsoluteTolerance;
    std::string mKspType;/**< KSP solver type (see PETSc KSPSetType() ) */
    std::string mPcType;/**< Preconditioner type (see PETSc PCSetType() ) */

    Vec mDirichletBoundaryConditionsVector; /**< Storage for efficient application of Dirichlet BCs, see AbstractBoundaryConditionsContainer */

    /** Stores a pointer to a purpose-build preconditioner*/
    PCBlockDiagonal* mpBlockDiagonalPC;

#ifdef TRACE_KSP
    unsigned mNumSolves;
    unsigned mTotalNumIterations;
    unsigned mMaxNumIterations;
#endif    
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        //MatNullSpace mMatNullSpace; ///\todo 
  
        archive & mNonZerosUsed;  
        archive & mMatrixIsConstant;
        archive & mTolerance; 
        archive & mUseAbsoluteTolerance;
        archive & mKspType;
        archive & mPcType;

        //Vec mDirichletBoundaryConditionsVector; ///\todo
    }    

public:

    /**
     * Constructor.
     *
     * @param lhsVectorSize
     * @param matType defaults to MATMPIAIJ
     */
    LinearSystem(PetscInt lhsVectorSize, MatType matType=(MatType) MATMPIAIJ);

    /**
     * Alternative constructor.
     *
     * Create a linear system, where the size is based on the size of a given
     * PETSc vec.
     *
     * The LHS & RHS vectors will be created by duplicating this vector's
     * settings.  This should avoid problems with using VecScatter on
     * bidomain simulation results.
     *
     * @param templateVector
     */
    LinearSystem(Vec templateVector);

    /**
     * Alternative constructor.
     *
     * Create a linear system which wraps the provided PETSc objects so we can
     * access them using our API.  Either of the objects may be NULL, but at
     * least one of them must not be.
     *
     * Useful for storing residuals and jacobians when solving nonlinear PDEs.
     *
     * @param residualVector
     * @param jacobianMatrix
     */
    LinearSystem(Vec residualVector, Mat jacobianMatrix);
    
    /**
     * Alternative constructor for archiving.
     * 
     * @param lhsVectorSize
     * @param lhsMatrix
     * @param rhsVector
     * @param matType defaults to MATMPIAIJ
     */
    LinearSystem(PetscInt lhsVectorSize, Mat lhsMatrix, Vec rhsVector, MatType matType=(MatType) MATMPIAIJ);

    /**
     * Destructor.
     */
    ~LinearSystem();
    
    /**
     * Helper method for the constructor. Initializes the LHS matrix and RHS vector.
     * 
     * @param matType
     */
    void SetupVectorAndMatrix(MatType matType);

//    bool IsMatrixEqualTo(Mat testMatrix);
//    bool IsRhsVectorEqualTo(Vec testVector);
    /**
     * Change one of the entires of the matrix to the specified value.
     * @param row
     * @param col
     * @param value
     */
    void SetMatrixElement(PetscInt row, PetscInt col, double value);
    /**
     * Add the specified value to an entry of the matrix.
     * @param row
     * @param col
     * @param value
     */
    void AddToMatrixElement(PetscInt row, PetscInt col, double value);

    /**
     * Call this before Solve().
     *
     * This calls AssembleFinalLhsMatrix() and AssembleRhsVector().
     */
    void AssembleFinalLinearSystem();
    /**
     * Should be called before AddToMatrixElement.
     *
     * This calls AssembleIntermediateLhsMatrix() and AssembleRhsVector().
     */
    void AssembleIntermediateLinearSystem();
    /**
     * Sets up the PETSc matrix left-hand-side mLhsMatrix
     */
    void AssembleFinalLhsMatrix();
    /**
     * Sets up the PETSc matrix left-hand-side mLhsMatrix
     */
    void AssembleIntermediateLhsMatrix();
    /**
     * Sets up the PETSc vector right-hand-side mRhsVector
     */
    void AssembleRhsVector();

    /**
     * Force PETSc to treat the matrix in this linear system as symmetric from now on.
     * @param isSymmetric Whether the matrix is symmetric or not
     */
    void SetMatrixIsSymmetric(bool isSymmetric=true);

    /**
     * Set mMatrixIsConstant.
     *
     * @param matrixIsConstant
     */
    void SetMatrixIsConstant(bool matrixIsConstant);

    /**
     * Set the relative tolerance.
     *
     * @param relativeTolerance
     */
    void SetRelativeTolerance(double relativeTolerance);

    /**
     * Set the absolute tolerance.
     *
     * @param absoluteTolerance
     */
    void SetAbsoluteTolerance(double absoluteTolerance);

    /**
     * Set the KSP solver type (see PETSc KSPSetType() for valid arguments)
     * @param kspType
     */
    void SetKspType(const char*);

    /**
     * Set the preconditioner type  (see PETSc PCSetType() for valid arguments)
     * @param pcType
     */
    void SetPcType(const char*);

    /**
     * Display the left-hand side matrix.
     */
    void DisplayMatrix();

    /**
     * Display the right-hand side vector.
     */
    void DisplayRhs();

    /**
     * Set all entries in a given row of a matrix to a certain value.
     * This must be called by the process who owns the row, (but other
     * processors will treat it as a null-op
     *
     * @param row
     * @param value
     */
    void SetMatrixRow(PetscInt row, double value);

    /**
     * Zero a row of the left-hand side matrix.
     * This method is a collective call (all processes should call it together).
     * If processes call it with different arguments then its results may 
     * not be predictable.
     * 
     * @param row
     */
    void ZeroMatrixRow(PetscInt row);

    /**
     * Zero a column of the left-hand side matrix.
     *
     * Unfortunately there is no equivalent method in Petsc, so this has to be
     * done carefully to ensure that the sparsity structure of the matrix
     * is not broken. Only owned entries which are non-zero are zeroed.
     *
     * @param col
     */
    void ZeroMatrixColumn(PetscInt col);

    /**
     * Zero all entries of the left-hand side matrix.
     */
    void ZeroLhsMatrix();

    /**
     * Zero all entries of the right-hand side vector.
     */
    void ZeroRhsVector();

    /**
     * Zero all entries of the left-hand side matrix and right-hand side vector.
     */
    void ZeroLinearSystem();

    /**
     * Solve the linear system.
     *
     * @param lhsGuess  an optional initial guess for the solution (defaults to NULL)
     */
    Vec Solve(Vec lhsGuess=NULL);

    /**
     * Set an element of the right-hand side vector to a given value.
     *
     * @param row
     * @param value
     */
    void SetRhsVectorElement(PetscInt row, double value);

    /**
     * Add a value to an element of the right-hand side vector.
     *
     * @param row
     * @param value
     */
    void AddToRhsVectorElement(PetscInt row, double value);

    /**
     * Get method for mSize.
     */
    unsigned GetSize() const;

    /**
     *
     * @param nullbasis
     * @param numberOfBases
     */
    void SetNullBasis(Vec nullbasis[], unsigned numberOfBases);

    /**
     * Get access to the rhs vector directly. Shouldn't generally need to be called.
     */
    Vec& rGetRhsVector();
    
    /**
     * Get access to the rhs vector for archiving
     */
    Vec GetRhsVector() const;
    
    /**
     * Get access to the lhs matrix directly. Shouldn't generally need to be called.
     */
    Mat& rGetLhsMatrix();
    
    /**
     * Get access to the lhs matrix for archiving
     */
    Mat GetLhsMatrix() const;

    /**
     * Gets access to the dirichlet boundary conditions vector.
     *
     * Should only be used by the BoundaryConditionsContainer.
     */
    Vec& rGetDirichletBoundaryConditionsVector();

    // DEBUGGING CODE:
    /**
     * Get this process's ownership range of the contents of the system.
     *
     * @param lo
     * @param hi
     */
    void GetOwnershipRange(PetscInt& lo, PetscInt& hi);

    /**
     * Return an element of the matrix.
     * May only be called for elements you own.
     *
     * @param row
     * @param col
     */
    double GetMatrixElement(PetscInt row, PetscInt col);

    /**
     * Return an element of the RHS vector.
     * May only be called for elements you own.
     *
     * @param row
     */
    double GetRhsVectorElement(PetscInt row);


    /**
     * Return the number of iterations taken by the last Solve()     
     */
    unsigned GetNumIterations() const;
    
    /**
     * Add multiple values to the matrix of linear system.
     *
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
        PetscInt num_rows_owned = 0;
        PetscInt global_row;

        for (unsigned row = 0; row<MATRIX_SIZE; row++)
        {
            global_row = matrixRowAndColIndices[row];
            if (global_row >=mOwnershipRangeLo && global_row <mOwnershipRangeHi)
            {
                matrix_row_indices[num_rows_owned++] = global_row;
            }
        }

        if ( num_rows_owned == MATRIX_SIZE)
        {
            MatSetValues(mLhsMatrix,
                         num_rows_owned,
                         matrix_row_indices,
                         MATRIX_SIZE,
                         (PetscInt*) matrixRowAndColIndices,
                         smallMatrix.data(),
                         ADD_VALUES);
        }
        else
        {
            // We need continuous data, if some of the rows do not belong to the processor their values
            // are not passed to MatSetValues
            double values[MATRIX_SIZE*MATRIX_SIZE];
            unsigned num_values_owned = 0;
            for (unsigned row = 0; row<MATRIX_SIZE; row++)
            {
                global_row = matrixRowAndColIndices[row];
                if (global_row >=mOwnershipRangeLo && global_row <mOwnershipRangeHi)
                {
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
        }
    };

    /**
     * Add multiple values to the RHS vector.
     *
     * @param vectorIndices mapping from index of the ublas vector (see param below)
     *  to index of the vector of this linear system
     * @param smallVector Ublas vector containing the values to be added
     *
     * N.B. Values which are not local (ie the row is not owned) will be skipped.
     */
    template<size_t VECTOR_SIZE>
    void AddRhsMultipleValues(unsigned* vectorIndices, c_vector<double, VECTOR_SIZE>& smallVector)
    {
        PetscInt indices_owned[VECTOR_SIZE];
        PetscInt num_indices_owned = 0;
        PetscInt global_row;

        for (unsigned row = 0; row<VECTOR_SIZE; row++)
        {
            global_row = vectorIndices[row];
            if (global_row >=mOwnershipRangeLo && global_row <mOwnershipRangeHi)
            {
                indices_owned[num_indices_owned++] = global_row;
            }
        }

        if (num_indices_owned == VECTOR_SIZE)
        {
            VecSetValues(mRhsVector,
                         num_indices_owned,
                         indices_owned,
                         smallVector.data(),
                         ADD_VALUES);
        }
        else
        {
            // We need continuous data, if some of the rows do not belong to the processor their values
            // are not passed to MatSetValues
            double values[VECTOR_SIZE];
            unsigned num_values_owned = 0;

            for (unsigned row = 0; row<VECTOR_SIZE; row++)
            {
                global_row = vectorIndices[row];
                if (global_row >=mOwnershipRangeLo && global_row <mOwnershipRangeHi)
                {
                    values[num_values_owned++] = smallVector(row);
                }
            }

            VecSetValues(mRhsVector,
                         num_indices_owned,
                         indices_owned,
                         values,
                         ADD_VALUES);
        }
    }

};

// Declare identifier for the serializer
BOOST_CLASS_EXPORT(LinearSystem);

namespace boost
{
namespace serialization
{

template<class Archive>
inline void save_construct_data(
    Archive & ar, const LinearSystem * t, const unsigned int file_version)
{
    
    OutputFileHandler handler("Archive", false);
    std::string archive_filename_lhs, archive_filename_rhs;
    archive_filename_lhs = handler.GetOutputDirectoryFullPath() + "lhs.arch";   
    archive_filename_rhs = handler.GetOutputDirectoryFullPath() + "rhs.arch"; 
    const unsigned size = t->GetSize();
    ar << size;
    
    PetscViewer vec_viewer;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, archive_filename_rhs.c_str(), FILE_MODE_WRITE, &vec_viewer);
    VecView(t->GetRhsVector(), vec_viewer);
    PetscViewerDestroy(vec_viewer);
    
    
    PetscViewer mat_viewer;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, archive_filename_lhs.c_str(), FILE_MODE_WRITE, &mat_viewer);
    MatView(t->GetLhsMatrix(), mat_viewer);
    PetscViewerDestroy(mat_viewer);

    //Is the matrix structurally symmetric?
    PetscTruth symm_set, is_symmetric;
    is_symmetric = PETSC_FALSE;
    //Note that the following call only changes is_symmetric when symm_set is true
    MatIsSymmetricKnown(t->GetLhsMatrix(), &symm_set, &is_symmetric);
    assert(symm_set == is_symmetric);
    ar << symm_set;
}    
    
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate an instance (using existing constructor)
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, LinearSystem * t, const unsigned int file_version)
{
     OutputFileHandler handler("Archive", false);
     std::string archive_filename_lhs, archive_filename_rhs;
     archive_filename_lhs = handler.GetOutputDirectoryFullPath() + "lhs.arch";   
     archive_filename_rhs = handler.GetOutputDirectoryFullPath() + "rhs.arch";
      
     PetscInt size; 
     ar >> size;
     
     PetscViewer vec_viewer;
     PetscViewerBinaryOpen(PETSC_COMM_WORLD, archive_filename_rhs.c_str(), FILE_MODE_READ, &vec_viewer);
     
     Vec new_vec;
     VecLoad(vec_viewer, PETSC_NULL, &new_vec);
     PetscViewerDestroy(vec_viewer);

     PetscViewer mat_viewer;
     PetscViewerBinaryOpen(PETSC_COMM_WORLD, archive_filename_lhs.c_str(), FILE_MODE_READ, &mat_viewer);
     Mat new_mat;

     MatLoad(mat_viewer, MATMPIAIJ, &new_mat);
     PetscViewerDestroy(mat_viewer);
     
     //This has to occur after the call to MatLoad as the matrix does not exist until MatLoad is called.
     //The property will be copied & set correctly in the LinearSystem constructor.
     PetscTruth symm_set;
     ar >> symm_set;
     if (symm_set == PETSC_TRUE)
     {
        MatSetOption(new_mat, MAT_SYMMETRIC);
        MatSetOption(new_mat, MAT_SYMMETRY_ETERNAL);
     }

     ::new(t)LinearSystem(size, new_mat, new_vec, MATMPIMAIJ);
}
}
} // namespace ...

#endif //_LINEARSYSTEM_HPP_
