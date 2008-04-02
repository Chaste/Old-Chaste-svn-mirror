/*
Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _ABSTRACTNONLINEARSOLVER_HPP_
#define _ABSTRACTNONLINEARSOLVER_HPP_

/**
 * Abstract Nonlinear equation system solver, dictates that each solver must have a
 * Solve function.
 */

#include <petscsnes.h>

class AbstractNonlinearSolver
{

public:
    virtual Vec Solve(PetscErrorCode (*pComputeResidual)(SNES,Vec,Vec,void*),
                      PetscErrorCode (*pComputeJacobian)(SNES,Vec,Mat*,Mat*,MatStructure*,void*),
                      Vec initialGuess,
                      void *pContext)=0;
                      
    virtual ~AbstractNonlinearSolver()
    {}
};

#endif // _ABSTRACTNONLINEARSOLVER_HPP_
