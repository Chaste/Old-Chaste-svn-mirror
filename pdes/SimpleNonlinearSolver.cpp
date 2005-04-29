/**
 * Concrete Simple Nonlinear PDE system solver.
 */

#include "SimpleNonlinearSolver.hpp"
#include <iostream>
#include "Exception.hpp"
#include <sstream>


/**
 * Simple Nonlinear PDE system solver, uses the PETSc SNES Solver, which uses Newton's method.
 * 
 * @param (*pComputeResidual)(SNES,Vec,Vec,void*) points to the function which 
 * computes the residual, it must take arguments SNES (a PETSc nonlinear solver 
 * object), Vec (current guess - a vector of the correct size), Vec (a Vec of the 
 * correct size in which the residual is returned), void* (a pointer to 
 * anything you may need to refer to when calculating the residual)
 * 
 * @param (*pComputeJacobian)(SNES,Vec,Mat*,Mat*,MatStructure*,void*) points to 
 * the function which computes the Jacobian, it must take arguments SNES (a PETSc 
 * nonlinear solver * object), Mat* (a pointer to the Jacobian matrix) ,Mat* (a pointer 
 * to a preconditioner matrix), MatStructure* (points to the PETSc matrix type e.g. AIJ), void* (a pointer to 
 * anything you may need to refer to when calculating the residual). 
 * 
 * @param residual A PETSc Vec of the correct size, elements do not need to be specified.
 * Used by SNES when calculating the residual at each iteration.
 * 
 * @param initialGuess A PETSc Vec of the correct size, containing initial guesses
 *  for the nonlinear solver.
 * 
 * @return Returns a PETSc Vec of the solution.
 * 
 * To be used in the form:
 * Vec answer;
 * Vec residual;
 * VecDuplicate(initialGuess,&residual);
 * VecDuplicate(initialGuess,&answer);
 * answer=solver->Solve(&ComputeResidual, &ComputeJacobian, residual, initialGuess, NULL);
 * 
 * In the same file, but outside this class the functions ComputeResidual and 
 * ComputeJacobian must sit, using the input arguments specified above.
 */
Vec SimpleNonlinearSolver::Solve(PetscErrorCode (*pComputeResidual)(SNES,Vec,Vec,void*),
                                 PetscErrorCode (*pComputeJacobian)(SNES,Vec,Mat*,Mat*,MatStructure*,void*),
                                 Vec residual, Vec initialGuess, void *pContext)
{
    SNES snes;
    Mat J; //Jacobian Matrix

    int N; //number of elements
    //get the size of the jacobian from the residual
    VecGetSize(residual,&N);
    MatCreate(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,N,N,&J);
    MatSetFromOptions(J);

    SNESCreate(PETSC_COMM_WORLD, &snes);
    SNESSetFunction(snes, residual, pComputeResidual, pContext);
    SNESSetJacobian(snes, J, J, pComputeJacobian, pContext);
    SNESSetType(snes,SNESLS);
    
    // x is the iteration vector SNES uses when solving, set equal to initialGuess to start with
    Vec x;
    VecDuplicate(initialGuess, &x);
    VecCopy(initialGuess, x);
    
    //std::cout << "Just about to call SNESSOlve" << std::endl << std::flush;
    SNESSolve(snes, x);
    SNESConvergedReason reason;
    SNESGetConvergedReason(snes,&reason);
    if (reason<0)
    {
    	std::stringstream reason_stream;
    	reason_stream << reason;
    	throw Exception("Nonlinear Solver did not converge. Petsc reason code:"
    	                +reason_stream.str()+" .");
    }

    return x;

}
