/**
 * Concrete Simple Nonlinear PDE system solver.
 */

#ifndef SIMPLENONLINEARSOLVER_HPP
#define SIMPLENONLINEARSOLVER_HPP

#include "petscsnes.h"
#include "AbstractNonlinearSolver.hpp"

class SimpleNonlinearSolver : public AbstractNonlinearSolver
{
	public:
    Vec Solve(PetscErrorCode (*pFunction)(SNES,Vec,Vec,void*),
              PetscErrorCode (*pComputeJacobian)(SNES,Vec,Mat*,Mat*,MatStructure*,void*),
              Vec residual, Vec initialGuess, void *context);

};
#endif
