#ifndef _LINEARSYSTEM_HPP_
#define _LINEARSYSTEM_HPP_

#include "petscvec.h"
#include "petscmat.h"
#include "AbstractLinearSolver.hpp"

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
	Vec mLhsVector;
	int mSize;
public:
    LinearSystem(int lhsVectorSize);
//    bool IsMatrixEqualTo(Mat testMatrix);
//    bool IsRhsVectorEqualTo(Vec testVector);
    void SetMatrixElement(int row, int col, double value);
    void AddToMatrixElement(int row, int col, double value);
    void AssembleFinalMatrix();         // Call before solve
    void AssembleIntermediateMatrix();  // Should be called before AddToMatrixElement
    void DisplayMatrix();
    void DisplayRhs() ;
    void SetMatrixRow(int row, double value);
    double GetMatrixElement(int row, int col);
    double GetRhsVectorElement(int index);
    Vec Solve(AbstractLinearSolver *pSolver);
    void SetRhsVectorElement(int row, double value); 
    void AddToRhsVectorElement(int row, double value);
	int GetSize();
};

#endif //_LINEARSYSTEM_HPP_
