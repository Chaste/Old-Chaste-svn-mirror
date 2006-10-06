#ifndef _SIMPLENONLINEARSOLVER_HPP_
#define _SIMPLENONLINEARSOLVER_HPP_

/**
 * Concrete Simple Nonlinear PDE system solver.
 */

#include <petscsnes.h>
#include "AbstractNonlinearSolver.hpp"

class SimplePetscNonlinearSolver : public AbstractNonlinearSolver
{
public:
    Vec Solve(PetscErrorCode (*pComputeResidual)(SNES,Vec,Vec,void*),
              PetscErrorCode (*pComputeJacobian)(SNES,Vec,Mat*,Mat*,MatStructure*,void*),
              Vec initialGuess, 
              void *pContext);
              
};

#endif // _SIMPLENONLINEARSOLVER_HPP_
