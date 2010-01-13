/*

Copyright (C) University of Oxford, 2005-2010

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


#include "SimpleNewtonNonlinearSolver.hpp"
#include <iostream>
#include <cassert>
#include "Exception.hpp"

SimpleNewtonNonlinearSolver::SimpleNewtonNonlinearSolver(double linearSolverRelativeTolerance)
        : mLinearSolverRelativeTolerance(linearSolverRelativeTolerance),
        mTolerance(1e-5),
        mWriteStats(false)
{
    mTestDampingValues.push_back(14);
    mTestDampingValues.push_back(-0.1);
    mTestDampingValues.push_back(0.05);
    for (unsigned i=1; i<=12; i++)
    {
        double val = double(i)/10;
        mTestDampingValues.push_back(val);
    }
}

SimpleNewtonNonlinearSolver::~SimpleNewtonNonlinearSolver()
{}

Vec SimpleNewtonNonlinearSolver::Solve(PetscErrorCode (*pComputeResidual)(SNES,Vec,Vec,void*),
                                       PetscErrorCode (*pComputeJacobian)(SNES,Vec,Mat*,Mat*,MatStructure*,void*),
                                       Vec initialGuess,
                                       void* pContext)
{
    PetscInt size;
    VecGetSize(initialGuess, &size);

    Vec current_solution;
    VecDuplicate(initialGuess, &current_solution);
    VecCopy(initialGuess, current_solution);

    LinearSystem linear_system(current_solution);

    (*pComputeResidual)(NULL, current_solution, linear_system.rGetRhsVector(), pContext);


    double residual_norm;
    VecNorm(linear_system.rGetRhsVector(), NORM_2, &residual_norm);
    double scaled_residual_norm = residual_norm/size;

    if (mWriteStats)
    {
        std::cout << "Newton's method:\n  Initial ||residual||/N = " << scaled_residual_norm
                  << "\n  Attempting to solve to tolerance " << mTolerance << "..\n";
    }

    double old_scaled_residual_norm;
    unsigned counter = 0;
    while (scaled_residual_norm > mTolerance)
    {
        counter++;

        // store the old norm to check with the new later
        old_scaled_residual_norm = scaled_residual_norm;

        // compute Jacobian and solve J dx = f for the (negative) update dx, (J the jacobian, f the residual)
        (*pComputeJacobian)(NULL, current_solution, &(linear_system.rGetLhsMatrix()), NULL, NULL, pContext);

        Vec negative_update = linear_system.Solve();


        Vec test_vec;
        VecDuplicate(initialGuess, &test_vec);

        double best_damping_factor = 1.0;
        double best_scaled_residual = 1e20; //large

        // loop over all the possible damping value and determine which gives
        // smallest residual
        for (unsigned i=0; i<mTestDampingValues.size(); i++)
        {
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
            // VecWAXPY(&a, x, y, w) computes w = ax + y
            PetscScalar alpha = -mTestDampingValues[i];
            VecWAXPY(&alpha, negative_update, current_solution, test_vec);
#else
            //[note: VecWAXPY(w,a,x,y) computes w = ax+y]
            VecWAXPY(test_vec,-mTestDampingValues[i],negative_update,current_solution);
#endif

            // compute new residual
            linear_system.ZeroLinearSystem();
            (*pComputeResidual)(NULL, test_vec, linear_system.rGetRhsVector(), pContext);
            VecNorm(linear_system.rGetRhsVector(), NORM_2, &residual_norm);
            scaled_residual_norm = residual_norm/size;

            if (scaled_residual_norm < best_scaled_residual)
            {
                best_scaled_residual = scaled_residual_norm;
                best_damping_factor = mTestDampingValues[i];
            }
        }
        VecDestroy(test_vec);

        // check the smallest residual was actually smaller than the
        // previous. If not, quit.
        if (best_scaled_residual > old_scaled_residual_norm)
        {
            // Free memory
            VecDestroy(current_solution);
            VecDestroy(negative_update);
            // Raise error
            std::stringstream error_message;
            error_message << "Iteration " << counter << ", unable to find damping factor such that residual decreases in update direction";
            EXCEPTION(error_message.str());
        }

        if (mWriteStats)
        {
            std::cout << "    Best damping factor = " << best_damping_factor << "\n";
        }


        // update solution: current_guess = current_solution - best_damping_factor*negative_update
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
        double minus_test_value = -best_damping_factor;
        VecAXPY(&minus_test_value, negative_update, current_solution);
#else
        //[note: VecAXPY(y,a,x) computes y = ax+y]
        VecAXPY(current_solution, -best_damping_factor, negative_update);
#endif
        scaled_residual_norm = best_scaled_residual;
        VecDestroy(negative_update);

        // compute best residual vector again and store in linear_system for next Solve()
        linear_system.ZeroLinearSystem();
        (*pComputeResidual)(NULL, current_solution, linear_system.rGetRhsVector(), pContext);

        if (mWriteStats)
        {
            std::cout << "    Iteration " << counter << ": ||residual||/N = " << scaled_residual_norm << "\n";
        }
    }

    if (mWriteStats)
    {
        std::cout << "  ..solved!\n\n";
    }

    return current_solution;
}


void SimpleNewtonNonlinearSolver::SetTolerance(double tolerance)
{
    assert(tolerance > 0);
    mTolerance = tolerance;
}

