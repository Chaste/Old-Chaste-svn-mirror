/*

Copyright (C) University of Oxford, 2008

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


#ifndef _TESTCARDIACUNEVEN_HPP_
#define _TESTCARDIACUNEVEN_HPP_


#include <cxxtest/TestSuite.h>
#include <petscvec.h>
#include <vector>
#include "PetscSetupAndFinalize.hpp"
#include "MonodomainProblem.hpp"
#include "BidomainProblem.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "CheckMonoLr91Vars.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "DistributedVector.hpp"
#include "PetscTools.hpp"
#include "EventHandler.hpp"
#include "HeartConfig.hpp"

class TestCardiacUneven : public CxxTest::TestSuite
{
public:
    void tearDown()
    {
        HeartConfig::Reset();   
    }

    // Solve on a 1D string of cells, 1mm long with a space step of 0.1mm.
    void TestMonodomainTwoUnevenProcessors()
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));        
        HeartConfig::Instance()->SetSimulationDuration(2); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("MonodomainUneven");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_1d");
                        
        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1> cell_factory;

        MonodomainProblem<1> monodomain_problem( &cell_factory );

        monodomain_problem.SetNodesPerProcessorFilename("heart/test/data/11_nodes_2_processors.txt");

        if(PetscTools::NumProcs() == 2)
        {
            monodomain_problem.Initialise();


            HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
            HeartConfig::Instance()->SetCapacitance(1.0);

            monodomain_problem.Solve();

            // test whether voltages and gating variables are in correct ranges
            CheckMonoLr91Vars<1>(monodomain_problem);

            // check some voltages
            ReplicatableVector voltage_replicated(monodomain_problem.GetVoltage());

            PetscInt petsc_lo, petsc_hi;
            VecGetOwnershipRange(monodomain_problem.GetVoltage(),&petsc_lo,&petsc_hi);

            if(PetscTools::GetMyRank() == 0)
            {
                TS_ASSERT_EQUALS(0, petsc_lo);
                TS_ASSERT_EQUALS(8, petsc_hi);
            }
            else
            {
                TS_ASSERT_EQUALS(8, petsc_lo);
                TS_ASSERT_EQUALS(11, petsc_hi);
            }

            double atol=5e-3;

            TS_ASSERT_DELTA(voltage_replicated[1], 20.7710232, atol);
            TS_ASSERT_DELTA(voltage_replicated[3], 21.5319692, atol);
            TS_ASSERT_DELTA(voltage_replicated[5], 22.9280817, atol);
            TS_ASSERT_DELTA(voltage_replicated[7], 24.0611303, atol);
            TS_ASSERT_DELTA(voltage_replicated[9], -0.770330519, atol);
            TS_ASSERT_DELTA(voltage_replicated[10], -19.2234919, atol);
        }
        else
        {
            TS_ASSERT_THROWS_ANYTHING(monodomain_problem.Initialise());
            EventHandler::Reset();
        }
    }
    // Solve on a 1D string of cells, 1mm long with a space step of 0.1mm.
    void TestBidomainTwoUnevenProcessors()
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));        
        HeartConfig::Instance()->SetSimulationDuration(2); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("BidomainUneven");
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_1d");
                        
        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1> cell_factory;

        BidomainProblem<1> bidomain_problem( &cell_factory );

        bidomain_problem.SetNodesPerProcessorFilename("heart/test/data/11_nodes_2_processors.txt");

        if(PetscTools::NumProcs() == 2)
        {
            bidomain_problem.Initialise();

            HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
            HeartConfig::Instance()->SetCapacitance(1.0);

            bidomain_problem.Solve();

            PetscInt petsc_lo, petsc_hi;
            VecGetOwnershipRange(bidomain_problem.GetVoltage(),&petsc_lo,&petsc_hi);

            if(PetscTools::GetMyRank() == 0)
            {
                TS_ASSERT_EQUALS(0, petsc_lo);
                TS_ASSERT_EQUALS(8*2, petsc_hi);
            }
            else
            {
                TS_ASSERT_EQUALS(8*2, petsc_lo);
                TS_ASSERT_EQUALS(11*2, petsc_hi);
            }
        }
        else
        {
            TS_ASSERT_THROWS_ANYTHING(bidomain_problem.Initialise());
            EventHandler::Reset();
        }
    }
};

#endif // _TESTCARDIACUNEVEN_HPP_
