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

#ifndef TESTCHECKPOINTING_HPP_
#define TESTCHECKPOINTING_HPP_

#include <cxxtest/TestSuite.h>
#include "PetscSetupAndFinalize.hpp"

#include "PlaneStimulusCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "BidomainProblem.hpp"

class TestCheckpointing : public CxxTest::TestSuite
{
public:

    /*
     *  This test makes sure that a simulation of x seconds can be run in multiple steps (SetSimulationTime() 
     * plus Solve()) and return the same results that a single call to Solve()  
     */
    void TestMultipleCallsToProblemSolve()
    {
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/3D_0_to_1mm_6000_elements");
        HeartConfig::Instance()->SetOutputDirectory("Monodomain3d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("monodomain3d");

        ///////////////////////////////////////////////////////////////////
        // Multiples calls to Solve()
        ///////////////////////////////////////////////////////////////////
        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 3> cell_factory(-600.0*1000);
        BidomainProblem<3> bidomain_problem_multiple( &cell_factory );

        bidomain_problem_multiple.Initialise();
        HeartConfig::Instance()->SetSimulationDuration(0.1);
        bidomain_problem_multiple.Solve();
        HeartConfig::Instance()->SetSimulationDuration(0.2);
        bidomain_problem_multiple.Solve();

        ///////////////////////////////////////////////////////////////////
        // Single call to Solve()
        ///////////////////////////////////////////////////////////////////
        HeartConfig::Instance()->SetOutputDirectory("Bidomain3d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain3d");

        BidomainProblem<3> bidomain_problem_single( &cell_factory );
        bidomain_problem_single.Initialise();
        HeartConfig::Instance()->SetSimulationDuration(0.2);        
        bidomain_problem_single.Solve();

        ///////////////////////////////////////////////////////////////////
        // compare
        ///////////////////////////////////////////////////////////////////
        DistributedVector single_solution = bidomain_problem_single.GetSolutionDistributedVector();
        DistributedVector multiple_solution = bidomain_problem_multiple.GetSolutionDistributedVector();
        DistributedVector::Stripe single_vm(single_solution,0);
        DistributedVector::Stripe single_phie(single_solution,1);
        DistributedVector::Stripe multiple_vm(multiple_solution,0);
        DistributedVector::Stripe multiple_phie(multiple_solution,1);
        for (DistributedVector::Iterator index = single_solution.Begin();
             index != single_solution.End();
             ++index)
        {
            TS_ASSERT_DELTA(single_vm[index], multiple_vm[index], 1e-8);
            TS_ASSERT_DELTA(single_phie[index], multiple_phie[index], 1e-8);
        }
        
    }
    
    
};

#endif /*TESTCHECKPOINTING_HPP_*/
