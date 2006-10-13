#ifndef SIMPLENEWTONNONLINEARSOLVER_HPP_
#define SIMPLENEWTONNONLINEARSOLVER_HPP_

#include "AbstractNonlinearSolver.hpp"
#include "SimpleLinearSolver.hpp"
#include "LinearSystem.hpp"

class SimpleNewtonNonlinearSolver : public AbstractNonlinearSolver
{
private :
    double mLinearSolverRelativeTolerance;
    double mTolerance;
    bool mWriteStats;

public :
    SimpleNewtonNonlinearSolver(double linearSolverRelativeTolerance = 1e-6);
    virtual ~SimpleNewtonNonlinearSolver();

    virtual Vec Solve(PetscErrorCode (*pComputeResidual)(SNES,Vec,Vec,void*),
                      PetscErrorCode (*pComputeJacobian)(SNES,Vec,Mat*,Mat*,MatStructure*,void*),
                      Vec initialGuess,
                      void *pContext);
    
    /*< Set a tolerance other than the default */
    void SetTolerance(double tolerance);
    
    void SetWriteStats(bool writeStats = true)
    {
        mWriteStats = writeStats;
    }
};

#endif /*SIMPLENEWTONNONLINEARSOLVER_HPP_*/
