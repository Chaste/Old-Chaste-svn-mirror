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

#ifndef TESTPCBLOCKDIAGONAL_HPP_
#define TESTPCBLOCKDIAGONAL_HPP_

#include <cxxtest/TestSuite.h>
#include "LinearSystem.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ReplicatableVector.hpp"
#include "Timer.hpp"


class TestPCBlockDiagonal : public CxxTest::TestSuite
{
public:

    void TestBasicFunctionality() throw (Exception)
    {
        unsigned system_size = 10;
        LinearSystem ls = LinearSystem(system_size);
        
        ls.SetAbsoluteTolerance(1e-6);
        ls.SetKspType("cg");
        ls.SetPcType("blockdiagonal");
        
        /// \todo: we want to load a mesh coming from a realistic problem instead
        for (unsigned row=0; row<system_size; row++)
        {
            ls.SetMatrixElement((PetscInt)row, (PetscInt)row, (double) row*1000 + row+1);

            for (unsigned col=0; col<row; col++)
            {
                ls.SetMatrixElement((PetscInt)row, (PetscInt)col, (double) row*10 + col+1);
                ls.SetMatrixElement((PetscInt)col, (PetscInt)row, (double) row*10 + col+1);
            }
        }
        ls.AssembleFinalLinearSystem();
        
        /// \todo: this rhs might be inconsistent (if it has components in the kernel of A). A problem if we solve the singular problem. Get one from a realistic simulation.              
        // Set b = A * [1 1 ... 1]'
        Vec all_ones = PetscTools::CreateVec(system_size, 1.0);
        MatMult(ls.GetLhsMatrix(), all_ones, ls.GetRhsVector());
        
        Vec solution = ls.Solve();
        
        ReplicatableVector solution_replicated;
        solution_replicated.ReplicatePetscVector(solution);
        
        for (unsigned row=0; row<system_size; row++)
        {
            TS_ASSERT_DELTA(solution_replicated[row], 1.0, 1e-6)
        }        
    }
    
    void TestBetterThanNoPreconditioning()
    {
        unsigned point_jacobi_its;
        unsigned block_diag_its;
        
        Timer::Reset();        
        {
            unsigned system_size = 100;
            LinearSystem ls = LinearSystem(system_size);
            
            ls.SetAbsoluteTolerance(1e-6);
            ls.SetKspType("cg");
            ls.SetPcType("none");
            
            Mat& system_matrix = ls.rGetLhsMatrix();
            PetscTools::ReadPetscObject(system_matrix, "notforrelease/test/data/matrices/cube_6000elems_half_activated.mat");
            
            Vec& system_rhs = ls.rGetRhsVector();
            PetscTools::ReadPetscObject(system_rhs, "notforrelease/test/data/matrices/cube_6000elems_half_activated.vec");
            
            ls.Solve();
            
            point_jacobi_its = ls.GetNumIterations();
        }        
        Timer::PrintAndReset("No preconditioning");
        
        {
            unsigned system_size = 100;
            LinearSystem ls = LinearSystem(system_size);
            
            ls.SetAbsoluteTolerance(1e-6);
            ls.SetKspType("cg");
            ls.SetPcType("blockdiagonal");

            Mat& system_matrix = ls.rGetLhsMatrix();
            PetscTools::ReadPetscObject(system_matrix, "notforrelease/test/data/matrices/cube_6000elems_half_activated.mat");
            
            Vec& system_rhs = ls.rGetRhsVector();
            PetscTools::ReadPetscObject(system_rhs, "notforrelease/test/data/matrices/cube_6000elems_half_activated.vec");
                        
            ls.Solve();
            
            block_diag_its = ls.GetNumIterations();
        }
        Timer::Print("Block diagonal preconditioner");

        std::cout << block_diag_its << " " << point_jacobi_its << std::endl;
        TS_ASSERT_LESS_THAN_EQUALS(block_diag_its, point_jacobi_its);
        
    }
};

#endif /*TESTPCBLOCKDIAGONAL_HPP_*/
