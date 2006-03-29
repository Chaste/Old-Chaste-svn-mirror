#ifndef _ABSTRACTNONLINEARSOLVER_HPP_
#define _ABSTRACTNONLINEARSOLVER_HPP_

/**
 * Abstract Nonlinear PDE system solver, dictates that each solver must have a
 * Solve function.
 */

#include <petscsnes.h>

class AbstractNonlinearSolver
{
    
public:
    virtual Vec Solve(PetscErrorCode (*pComputeResidual)(SNES,Vec,Vec,void*), 
                      PetscErrorCode (*pComputeJacobian)(SNES,Vec,Mat*,Mat*,MatStructure*,void*),
                      Vec residual, Vec initialGuess,
                      void *context)=0;
    virtual ~AbstractNonlinearSolver()
    {
    }                  
};

#endif // _ABSTRACTNONLINEARSOLVER_HPP_
