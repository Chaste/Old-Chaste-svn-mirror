#ifndef _SIMPLELINEARSOLVER_H_
#define _SIMPLELINEARSOLVER_H_

#include "AbstractLinearSolver.hpp"
#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"

class SimpleLinearSolver : public AbstractLinearSolver
{

public:

    Vec Solve(Mat lhsMatrix, Vec rhsVector, int size);
    SimpleLinearSolver()
    {
        mLinearSystemKnown=false;
        mMatrixIsConstant=false;
    }
    ~SimpleLinearSolver()
    {
         if (mLinearSystemKnown==true)
         {
            KSPDestroy(mSimpleSolver);
         }
    }
    
    void SetMatrixIsConstant(){
        mMatrixIsConstant=true;
    }    
private:
    bool mLinearSystemKnown;
    bool mMatrixIsConstant;
    KSP mSimpleSolver;
    
};

#endif // _SIMPLELINEARSOLVER_H_
