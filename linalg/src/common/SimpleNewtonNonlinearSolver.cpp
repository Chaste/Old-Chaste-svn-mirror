#include <iostream>
#include "SimpleNewtonNonlinearSolver.hpp"
#include "Exception.hpp"

SimpleNewtonNonlinearSolver::SimpleNewtonNonlinearSolver()
{
    mpLinearSolver = new SimpleLinearSolver;
    mWeAllocatedSolverMemory = true;
    mTolerance = 1e-5;
    mWriteStats = false;
}

SimpleNewtonNonlinearSolver::~SimpleNewtonNonlinearSolver()
{
    if(mWeAllocatedSolverMemory)
    {
        delete mpLinearSolver;
    }
}

Vec SimpleNewtonNonlinearSolver::Solve(PetscErrorCode (*pComputeResidual)(SNES,Vec,Vec,void*),
                                       PetscErrorCode (*pComputeJacobian)(SNES,Vec,Mat*,Mat*,MatStructure*,void*),
                                       Vec initialGuess,
                                       void* pContext)
{
    int size;
    VecGetSize(initialGuess, &size);
    
    Vec current_solution;     
    VecDuplicate(initialGuess, &current_solution);
    VecCopy(initialGuess, current_solution);

    // don't need initialGuess anymore
    VecDestroy(initialGuess);
    
    LinearSystem linear_system(current_solution);
    
    (*pComputeResidual)(NULL, current_solution, linear_system.rGetRhsVector(), pContext);
    
    
    double residual_norm;
    VecNorm(linear_system.rGetRhsVector(), NORM_2, &residual_norm);
    double scaled_norm = residual_norm/size;
    
    if(mWriteStats)
    {
        std::cout << "Newton's method:\n  Initial ||residual||/N = " << scaled_norm 
                  << "\n  Attempting to solve to tolerance " << mTolerance << "..\n";
    }


    double old_scaled_norm;
    
    unsigned counter = 0;
    while(scaled_norm > mTolerance)
    {
        counter++;
        
        // store the old norm to check with the new later
        old_scaled_norm = scaled_norm;

        // compute Jacobian and solve J dx = f for the (negative) update dx, (J the jacobian, f the residual) 
        (*pComputeJacobian)(NULL, current_solution, &(linear_system.rGetLhsMatrix()), NULL, NULL, pContext);
     
        Vec update = linear_system.Solve(mpLinearSolver);

        // update to solution: current_guess += -update [note: VecAXPY(y,a,x) computes y = ax+y]
        VecAXPY(current_solution, -1, update);

        // compute new residual and check it didn't decrease
        linear_system.ZeroLinearSystem();
        (*pComputeResidual)(NULL, current_solution, linear_system.rGetRhsVector(), pContext);
        VecNorm(linear_system.rGetRhsVector(), NORM_2, &residual_norm);
        scaled_norm = residual_norm/size;
    
        if(scaled_norm > old_scaled_norm)
        {
            std::stringstream error_message;
            error_message << "Residual norm increased on iteration " << counter;
            EXCEPTION(error_message.str());
        }
        
        if(mWriteStats)
        {
            std::cout << "    Iteration " << counter << ": ||residual||/N = " << scaled_norm << "\n";
        }
    }


    if(mWriteStats)
    {
        std::cout << "  ..solved!\n\n";
    }

    return current_solution;
}



void SimpleNewtonNonlinearSolver::SetLinearSolver(AbstractLinearSolver* pLinearSolver)
{
    if(mWeAllocatedSolverMemory)
    {
        delete mpLinearSolver;
    }
    mpLinearSolver = pLinearSolver;
    mWeAllocatedSolverMemory = false;
}

void SimpleNewtonNonlinearSolver::SetTolerance(double tolerance)
{
    assert(tolerance > 0);
    mTolerance = tolerance;
}
    
