#ifndef _SIMPLENONLINEARSOLVER_HPP_
#define _SIMPLENONLINEARSOLVER_HPP_

/**
 * Concrete Simple Nonlinear PDE system solver.
 */

#include <petscsnes.h>
#include "AbstractNonlinearSolver.hpp"

class SimpleNonlinearSolver : public AbstractNonlinearSolver
{
public:
    Vec Solve(PetscErrorCode (*pFunction)(SNES,Vec,Vec,void*),
              PetscErrorCode (*pComputeJacobian)(SNES,Vec,Mat*,Mat*,MatStructure*,void*),
              Vec residual, Vec initialGuess, void *context);
              
};

#endif // _SIMPLENONLINEARSOLVER_HPP_
