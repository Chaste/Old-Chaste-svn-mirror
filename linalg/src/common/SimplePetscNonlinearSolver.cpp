/*

Copyright (C) University of Oxford, 2005-2009

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


/**
 * Concrete Simple Nonlinear PDE system solver.
 */



#include "SimplePetscNonlinearSolver.hpp"
#include "Exception.hpp"
#include "petscsnes.h"
#include "PetscTools.hpp"
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
                                      void* pContext)
{
    SNES snes;

    // create the residual vector by copying the structure of the initial guess
    Vec residual;
    VecDuplicate(initialGuess, &residual);

    Mat jacobian; //Jacobian Matrix

    PetscInt N; //number of elements
    //get the size of the jacobian from the residual
    VecGetSize(initialGuess,&N);

    PetscTools::SetupMat(jacobian, N, N);

    SNESCreate(PETSC_COMM_WORLD, &snes);
    SNESSetFunction(snes, residual, pComputeResidual, pContext);
    SNESSetJacobian(snes, jacobian, jacobian, pComputeJacobian, pContext);
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
    MatDestroy(jacobian); // Free Jacobian

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
