/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef TESTPCBLOCKDIAGONAL_HPP_
#define TESTPCBLOCKDIAGONAL_HPP_

#include <cxxtest/TestSuite.h>
#include "LinearSystem.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ReplicatableVector.hpp"
#include "Timer.hpp"
#include "DistributedVectorFactory.hpp"
#include <cstring>


/*
 *  Warning: these tests do not inform PETSc about the nullspace of the matrix. Therefore, convergence might be
 * different compared to a real cardiac simulation. Do not take conclusions about preconditioner/solver performance
 * based on these tests only. 
 */
class TestChebyshevIteration : public CxxTest::TestSuite
{
public:

    void TestChebyshevVsCG() throw (Exception)
    {
        unsigned num_nodes = 1331;
        DistributedVectorFactory factory(num_nodes);
        Vec parallel_layout = factory.CreateVec(2);

        unsigned cg_its;
        unsigned chebyshev_its;

        Timer::Reset();
        {
            Mat system_matrix;
            //Note that this test deadlocks if the file's not on the disk
            PetscTools::ReadPetscObject(system_matrix, "linalg/test/data/matrices/cube_6000elems_half_activated.mat", parallel_layout);

            Vec system_rhs;
            //Note that this test deadlocks if the file's not on the disk
            PetscTools::ReadPetscObject(system_rhs, "linalg/test/data/matrices/cube_6000elems_half_activated.vec", parallel_layout);

            LinearSystem ls = LinearSystem(system_rhs, system_matrix);

            ls.SetMatrixIsSymmetric();
            ls.SetAbsoluteTolerance(1e-9);
            ls.SetKspType("cg");
            ls.SetPcType("bjacobi");

            Vec solution = ls.Solve();

            cg_its = ls.GetNumIterations();

            MatDestroy(system_matrix);
            VecDestroy(system_rhs);
            VecDestroy(solution);
        }
        Timer::PrintAndReset("CG");

        {
            Mat system_matrix;
            //Note that this test deadlocks if the file's not on the disk
            PetscTools::ReadPetscObject(system_matrix, "linalg/test/data/matrices/cube_6000elems_half_activated.mat", parallel_layout);

            Vec system_rhs;
            //Note that this test deadlocks if the file's not on the disk
            PetscTools::ReadPetscObject(system_rhs, "linalg/test/data/matrices/cube_6000elems_half_activated.vec", parallel_layout);

            LinearSystem ls = LinearSystem(system_rhs, system_matrix);

            ls.SetMatrixIsSymmetric();
            ls.SetAbsoluteTolerance(1e-9);
            ls.SetKspType("chebychev");
            ls.SetPcType("bjacobi");

            Vec solution = ls.Solve();

            chebyshev_its = ls.GetNumIterations();

            MatDestroy(system_matrix);
            VecDestroy(system_rhs);
            VecDestroy(solution);
        }
        Timer::Print("Chebyshev");

        TS_ASSERT_LESS_THAN(cg_its, 15u); // Takes 14 iterations with 16 cores
        TS_ASSERT_LESS_THAN(chebyshev_its, 17u); // Takes 16 iterations with 16 cores

        VecDestroy(parallel_layout);
    }

    void TestChebyshevNoSpectrumShift() throw (Exception)
    {
        unsigned num_nodes = 1331;
        DistributedVectorFactory factory(num_nodes);
        Vec parallel_layout = factory.CreateVec(2);

        Vec zero_guess = factory.CreateVec(2);
        double zero = 0.0;
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2)
        VecSet(&zero, zero_guess);
#else
        VecSet(zero_guess, zero);
#endif
        unsigned chebyshev_its;

        Timer::Reset();

        Mat system_matrix;
        //Note that this test deadlocks if the file's not on the disk
        PetscTools::ReadPetscObject(system_matrix, "linalg/test/data/matrices/cube_6000elems_half_activated.mat", parallel_layout);

        Vec system_rhs;
        //Note that this test deadlocks if the file's not on the disk
        PetscTools::ReadPetscObject(system_rhs, "linalg/test/data/matrices/cube_6000elems_half_activated.vec", parallel_layout);

        LinearSystem ls = LinearSystem(system_rhs, system_matrix);

        ls.SetMatrixIsSymmetric();

        // Solve to relative convergence for coverage
        ls.SetRelativeTolerance(1e-9);
        ls.SetKspType("chebychev");

        // Solve with Jacobi to cover not using spectrum shift
        ls.SetPcType("jacobi");

        // Solving with zero guess for coverage.
        Vec solution = ls.Solve(zero_guess);

        chebyshev_its = ls.GetNumIterations();

        MatDestroy(system_matrix);
        VecDestroy(system_rhs);
        VecDestroy(solution);
        Timer::PrintAndReset("Chebyshev-Jacobi");

        // The number of iterations is PETSc-version-dependent...
        TS_ASSERT_LESS_THAN(130u, chebyshev_its);
        TS_ASSERT_LESS_THAN(chebyshev_its, 150u);

        VecDestroy(parallel_layout);
        VecDestroy(zero_guess);
    }

};

#endif /*TESTPCBLOCKDIAGONAL_HPP_*/
