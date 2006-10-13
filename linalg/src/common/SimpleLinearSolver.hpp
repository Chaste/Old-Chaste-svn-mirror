#ifndef _SIMPLELINEARSOLVER_H_
#define _SIMPLELINEARSOLVER_H_

#include "AbstractLinearSolver.hpp"
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>

class SimpleLinearSolver : public AbstractLinearSolver
{
private:
    double mRelativeTolerance;

public:

    Vec Solve(Mat lhsMatrix, Vec rhsVector, int size, MatNullSpace matNullSpace=NULL);
    SimpleLinearSolver(double relTolerance): mRelativeTolerance(relTolerance)
    {
        mLinearSystemKnown=false;
        mMatrixIsConstant=false;
    }
    virtual ~SimpleLinearSolver()
    {
        if (mLinearSystemKnown==true)
        {
            KSPDestroy(mSimpleSolver);
        }
    }
    
    void SetMatrixIsConstant()
    {
        mMatrixIsConstant=true;
    }
    
private:
    bool mLinearSystemKnown;
    bool mMatrixIsConstant;
    KSP mSimpleSolver;
    
};

#endif // _SIMPLELINEARSOLVER_H_
