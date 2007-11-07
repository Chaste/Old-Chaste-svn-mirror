#ifndef SIMPLENEWTONNONLINEARSOLVER_HPP_
#define SIMPLENEWTONNONLINEARSOLVER_HPP_

#include "LinearSystem.hpp"
#include "AbstractNonlinearSolver.hpp"
#include "SimpleLinearSolver.hpp"
#include <vector>

class SimpleNewtonNonlinearSolver : public AbstractNonlinearSolver
{
private :
    double mLinearSolverRelativeTolerance;
    double mTolerance;
    bool mWriteStats;
    
    std::vector<double> mTestDampingValues;
    
public :
    SimpleNewtonNonlinearSolver(double linearSolverRelativeTolerance = 1e-6);
    virtual ~SimpleNewtonNonlinearSolver();
    
    /**
     * Solve()
     * 
     * Solve a nonlinear system using Newton's method with damping. Newton's algorithm
     * is 
     * 
     * x_new = x_old - J^{-1} f
     * 
     * where J is the Jacobian matrix evaluated at x_old and f the residual evaluated at 
     * x_old. The Newton method with damping is 
     * 
     * x_new = x_old - s J^{-1} f
     * 
     * where s is some damping factor. Here s is chosen by just looked at a fixed set of 
     * possible damping factors and choosing the one which gives the best x_new (the one 
     * for which the residual evaluated at x_new has the lowest norm).
     * 
     * The solver quits once the ||f||/numVariables 
     * 
     * @param pComputeResidual points to the function which
     * computes the residual, it must take arguments SNES (a PETSc nonlinear solver
     * object), Vec (current guess - a vector of the correct size), Vec (a Vec of the
     * correct size in which the residual is returned), void* (a pointer to
     * anything you may need to refer to when calculating the residual)
     *
     * @param pComputeJacobian points to
     * the function which computes the Jacobian, it must take arguments SNES (a PETSc
     * nonlinear solver * object), Mat* (a pointer to the Jacobian matrix) ,Mat* (a pointer
     * to a preconditioner matrix), MatStructure* (points to the PETSc matrix type e.g. AIJ), void* (a pointer to
     * anything you may need to refer to when calculating the residual).
     * 
     * @param initialGuess A PETSc Vec of the correct size, containing initial guesses
     * for the nonlinear solver.
     * 
     * @param pContext [optional] A pointer to a class that may have to be used in the 
     *  ComputeResidual and ComputeJacobian functions
     *
     * @return Returns a PETSc Vec of the solution.
     *
     * To be used in the form:
     * Vec answer=solver->Solve(&ComputeResidual, &ComputeJacobian, initialGuess, NULL);
     *
     * In the same file, but outside this class the functions ComputeResidual and
     * ComputeJacobian must sit, using the input arguments specified above.
     */
    virtual Vec Solve(PetscErrorCode (*pComputeResidual)(SNES,Vec,Vec,void*),
                      PetscErrorCode (*pComputeJacobian)(SNES,Vec,Mat*,Mat*,MatStructure*,void*),
                      Vec initialGuess,
                      void *pContext);
                      
    /*< Set a tolerance other than the default */
    void SetTolerance(double tolerance);
    
    /*< Call to set the solver to write details as it solves */
    void SetWriteStats(bool writeStats = true)
    {
        mWriteStats = writeStats;
    }
};

#endif /*SIMPLENEWTONNONLINEARSOLVER_HPP_*/
