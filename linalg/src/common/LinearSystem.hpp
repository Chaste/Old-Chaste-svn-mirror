#ifndef _LINEARSYSTEM_HPP_
#define _LINEARSYSTEM_HPP_

#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include "UblasCustomFunctions.hpp"

class AbstractLinearSolver;

#include <string>

/**
 * Linear System class. Stores and solves a linear equation of the form Ax=b,
 * where A is a square matrix and x and b are column vectors.
 * The class uses PETSc.
 */
class LinearSystem
{
private:
    Mat mLhsMatrix;
    Vec mRhsVector;
    //Vec mLhsVector;
    PetscInt mSize;
    /** \todo
     * Verify claim that ownership range for Vec and Mat is same.
     * This should only matter for efficiency if the claim is false.
     */
    PetscInt mOwnershipRangeLo; /*< For parallel code.  Stores lowest index of vectors and lowest row of matrix*/
    PetscInt mOwnershipRangeHi; /*< Stores <b>one more than</b> the highest index stored locally*/
    
    MatNullSpace mMatNullSpace;
    
    /** Whether we need to destroy the PETSc objects in our destructor */
    bool mDestroyPetscObjects;
    
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
    
    void DisplayMatrix();
    void DisplayRhs() ;
    void SetMatrixRow(PetscInt row, double value);
    void ZeroMatrixRow(PetscInt row);
    void ZeroLhsMatrix();
    void ZeroRhsVector();
    void ZeroLinearSystem();
    Vec Solve(AbstractLinearSolver *pSolver, Vec lhsGuess=NULL);
    void SetRhsVectorElement(PetscInt row, double value);
    void AddToRhsVectorElement(PetscInt row, double value);
    void AddToRhsVectorElements(PetscInt m_rhs, PetscInt idx_rhs[], double rhs[]);
    unsigned GetSize();
    void SetNullBasis(Vec nullbasis[], unsigned numberOfBases);
    Vec& rGetRhsVector();
    Mat& rGetLhsMatrix();
    
    
    // DEBUGGING CODE:
    void GetOwnershipRange(PetscInt &lo, PetscInt &hi);
    double GetMatrixElement(PetscInt row, PetscInt col);
    double GetRhsVectorElement(PetscInt row);
    //void WriteLinearSystem(std::string matrixFile, std::string rhsVectorFile);
    
    
    /***
     * Add multiple values to a the linear system
     * @param matrixRowAndColIndices mapping from index of the ublas matrix (see param below)
     *  to index of the Petsc matrix of this linear system
     * @param smallMatrix Ublas matrix containing the values to be added
     * 
     * N.B. Values which are not local (ie the row is not owned) will be skipped.
     */
      
    template<unsigned MAT_SIZE>
    void AddMultipleValues(unsigned* matrixRowAndColIndices, c_matrix<double, MAT_SIZE, MAT_SIZE>& smallMatrix)
    {
        PetscInt matrix_row_indices[MAT_SIZE];
        PetscInt num_rows_owned=0; 
        unsigned num_values_owned=0;
        
        double values[MAT_SIZE*MAT_SIZE];
        for (unsigned row = 0 ; row<MAT_SIZE; row++)
        {
            PetscInt global_row = matrixRowAndColIndices[row];
            if (global_row >=mOwnershipRangeLo && global_row <mOwnershipRangeHi)
            {
                matrix_row_indices[num_rows_owned++] = global_row ;
                for (unsigned col=0; col<MAT_SIZE; col++)
                {
                    values[num_values_owned++] = smallMatrix(row,col);
                }
            }
        }
        
        
    
        MatSetValues(mLhsMatrix,
                     num_rows_owned,
                     matrix_row_indices,
                     MAT_SIZE,
                     (PetscInt*) matrixRowAndColIndices,
                     values,
                     ADD_VALUES);
        
    };
    
    
};

#endif //_LINEARSYSTEM_HPP_
