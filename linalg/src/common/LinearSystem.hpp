#ifndef _LINEARSYSTEM_HPP_
#define _LINEARSYSTEM_HPP_

#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>


class AbstractLinearSolver;


#include <string>

/**
 * Linear System class. Stores and solves a linear equation of the form Ax=b,
 * where A is a square matrix and x and b are column vectors.
 * The class uses PETSc.
 *
 *
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
    
public:
    LinearSystem(PetscInt lhsVectorSize);
    LinearSystem(Vec templateVector);
    ~LinearSystem();
//    bool IsMatrixEqualTo(Mat testMatrix);
//    bool IsRhsVectorEqualTo(Vec testVector);
    void SetMatrixElement(PetscInt row, PetscInt col, double value);
    void AddToMatrixElement(PetscInt row, PetscInt col, double value);
    void AssembleFinalLinearSystem();         // Call before solve
    void AssembleIntermediateLinearSystem();  // Should be called before AddToMatrixElement
    void AssembleRhsVector();
    void DisplayMatrix();
    void DisplayRhs() ;
    void SetMatrixRow(PetscInt row, double value);
    void ZeroMatrixRow(PetscInt row);
    void ZeroLhsMatrix();
    void ZeroRhsVector();
    void ZeroLinearSystem();
    Vec Solve(AbstractLinearSolver *pSolver);
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
    //void WriteLinearSystem(std::string matrixFile, std::string rhsVectorFile);
    
};

#endif //_LINEARSYSTEM_HPP_
