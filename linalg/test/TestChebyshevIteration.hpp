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
 * different compared to a real cardiac simulation. Do not take conclusions about preconditioner performance
 * based on these tests only. 
 */
class TestChebyshevIteration : public CxxTest::TestSuite
{
public:

    void TestChebyshevVsCG()
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
            ls.SetPcType("jacobi");

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
            ls.SetPcType("jacobi");

            Vec solution = ls.Solve();

            chebyshev_its = ls.GetNumIterations();

            MatDestroy(system_matrix);
            VecDestroy(system_rhs);
            VecDestroy(solution);
        }
        Timer::Print("Chebyshev");

        TS_ASSERT_EQUALS(cg_its, 101u);
        TS_ASSERT_EQUALS(chebyshev_its, 175u);

        VecDestroy(parallel_layout);
    }

    void TestFixedNumberOfIterations()
    {
        unsigned num_nodes = 1331;
        DistributedVectorFactory factory(num_nodes);
        Vec parallel_layout = factory.CreateVec(2);               

        Mat system_matrix;
        //Note that this test deadlocks if the file's not on the disk
        PetscTools::ReadPetscObject(system_matrix, "linalg/test/data/matrices/cube_6000elems_half_activated.mat", parallel_layout);

        Vec system_rhs;
        //Note that this test deadlocks if the file's not on the disk
        PetscTools::ReadPetscObject(system_rhs, "linalg/test/data/matrices/cube_6000elems_half_activated.vec", parallel_layout);

        LinearSystem ls = LinearSystem(system_rhs, system_matrix);

        ls.SetMatrixIsSymmetric();
        ls.SetKspType("chebychev");
        ls.SetPcType("jacobi");
        ls.SetUseFixedNumberIterations();

        Vec guess;
        VecDuplicate(system_rhs, &guess);
        VecSet(guess, 0.0);

        Vec solution = ls.Solve(guess);

        unsigned chebyshev_its = ls.GetNumIterations();
        TS_ASSERT_EQUALS(chebyshev_its, 88u);

        /*
         * We solve the same linear system again using the previous solution as the new
         * guess. If we were checking convergence normally it would take one iteration
         * to solve. Since we set fixed number of iterations based on first solve, it will
         * take the same number as above.
         */
        Vec solution2 = ls.Solve(solution);
        chebyshev_its = ls.GetNumIterations();
        TS_ASSERT_EQUALS(chebyshev_its, 88u);

        MatDestroy(system_matrix);
        VecDestroy(system_rhs);
        VecDestroy(solution);
        VecDestroy(solution2);
        
        VecDestroy(parallel_layout);
    }

};

#endif /*TESTPCBLOCKDIAGONAL_HPP_*/
