/**
 * Concrete Simple Nonlinear PDE system solver.
 */



#include "SimplePetscNonlinearSolver.hpp"
#include "Exception.hpp"
#include "petscsnes.h"
#include <sstream>


/**
 * Simple Nonlinear PDE system solver, uses the PETSc SNES Solver, which uses Newton's method.
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
 *  for the nonlinear solver.
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
Vec SimplePetscNonlinearSolver::Solve(PetscErrorCode (*pComputeResidual)(SNES,Vec,Vec,void*),
                                 PetscErrorCode (*pComputeJacobian)(SNES,Vec,Mat*,Mat*,MatStructure*,void*),
                                 Vec initialGuess, 
                                 void *pContext)
{
    SNES snes;
    
    // create the residual vector by copying the structure of the initial guess
    Vec residual;
    VecDuplicate(initialGuess, &residual);
    
    Mat J; //Jacobian Matrix
    
    int N; //number of elements
    //get the size of the jacobian from the residual
    VecGetSize(initialGuess,&N);
    
    
#if (PETSC_VERSION_MINOR == 2) //Old API
    MatCreate(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,N,N,&J);
#else
    MatCreate(PETSC_COMM_WORLD,&J);
    MatSetSizes(J,PETSC_DECIDE,PETSC_DECIDE,N,N);
#endif
    MatSetFromOptions(J);
    
    SNESCreate(PETSC_COMM_WORLD, &snes);
    SNESSetFunction(snes, residual, pComputeResidual, pContext);
    SNESSetJacobian(snes, J, J, pComputeJacobian, pContext);
    SNESSetType(snes,SNESLS);
    SNESSetTolerances(snes,1.0e-5,1.0e-5,1.0e-5,PETSC_DEFAULT,PETSC_DEFAULT);
    
    // x is the iteration vector SNES uses when solving, set equal to initialGuess to start with
    Vec x;
    VecDuplicate(initialGuess, &x);
    VecCopy(initialGuess, x);
    
    
#if (PETSC_VERSION_MINOR == 2) //Old API
    SNESSolve(snes, x);
#else
    SNESSolve(snes, PETSC_NULL, x);
#endif
    
    VecDestroy(residual);
    MatDestroy(J); // Free Jacobian
    
    SNESConvergedReason reason;
    SNESGetConvergedReason(snes,&reason);
#define COVERAGE_IGNORE
    if (reason<0)
    {
        std::stringstream reason_stream;
        reason_stream << reason;
        VecDestroy(x); // Since caller can't free the memory in this case
        SNESDestroy(snes);
        EXCEPTION("Nonlinear Solver did not converge. Petsc reason code:"
                  +reason_stream.str()+" .");
    }
#undef COVERAGE_IGNORE
    SNESDestroy(snes);
    
    return x;
}
