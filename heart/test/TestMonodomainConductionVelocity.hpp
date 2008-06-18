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


#ifndef _TESTMONODOMAINCONDUCTIONVELOCITY_HPP_
#define _TESTMONODOMAINCONDUCTIONVELOCITY_HPP_


#include <cxxtest/TestSuite.h>
#include <vector>

#include "ConformingTetrahedralMesh.hpp"
#include "PropagationPropertiesCalculator.hpp"
#include "Hdf5DataReader.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "MonodomainProblem.hpp"
#include "CheckMonoLr91Vars.hpp"
#include "PlaneStimulusCellFactory.hpp"


class TestMonodomainConductionVelocity : public CxxTest::TestSuite
{
public:
    void tearDown()
    {
        HeartConfig::Destroy();   
    }

    // Solve on a 1D string of cells, 1cm long with a space step of 0.1mm.
    void TestMonodomainDg01DWith100elements()
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        
        PlaneStimulusCellFactory<1> cell_factory;
        MonodomainProblem<1> monodomain_problem(&cell_factory);

        monodomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1_100_elements");
        monodomain_problem.SetEndTime(30);   // 30 ms
        monodomain_problem.SetOutputDirectory("MonoConductionVel");
        monodomain_problem.SetOutputFilenamePrefix("MonodomainLR91_1d");

        std::vector<unsigned> output_nodes;
        output_nodes.push_back(5);
        output_nodes.push_back(95);
        monodomain_problem.SetOutputNodes(output_nodes);

        monodomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        monodomain_problem.Solve();

        // test whether voltages and gating variables are in correct ranges
        CheckMonoLr91Vars<1>(monodomain_problem);

        // Calculate the conduction velocity
        Hdf5DataReader simulation_data("MonoConductionVel",
                                       "MonodomainLR91_1d");
        PropagationPropertiesCalculator ppc(&simulation_data);
        double velocity=0.0;

        // Check action potential propagated to node 95
        TS_ASSERT_THROWS_NOTHING(velocity=ppc.CalculateConductionVelocity(5,95,0.9));

        // The value should be approximately 50cm/sec
        // i.e. 0.05 cm/msec (which is the units of the simulation)
        TS_ASSERT_DELTA(velocity, 0.05, 0.003);
    }

    // Solve on a 1D string of cells, 1cm long with a space step of 0.5mm.
    //
    // Note that this space step ought to be too big!
    void TestMonodomainDg01DWith20elements()
    {
//#ifdef NDEBUG
//        TS_ASSERT(true);
//#else
        //Note that this test *used to* FAIL under any NDEBUG build.
        /*
         * This is because it is testing exceptions which are tripped by gating variables going
         * out of range in the Luo-Rudy cell model.  These variable ranges are not tested in
         * production builds.  They are guarded with  "#ifndef NDEBUG"
         */

        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));

        PlaneStimulusCellFactory<1> cell_factory;
        MonodomainProblem<1> monodomain_problem(&cell_factory);

        monodomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1_20_elements");
        monodomain_problem.SetEndTime(1);   // 1 ms
        monodomain_problem.SetOutputDirectory("MonoConductionVel");
        monodomain_problem.SetOutputFilenamePrefix("MonodomainLR91_1d");

        monodomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        // the mesh is too coarse, and this simulation will result in cell gating
        // variables going out of range. An exception should be thrown in the
        // EvaluateYDerivatives() method of the cell model

        TS_ASSERT_THROWS_ANYTHING(monodomain_problem.Solve());
//#endif
    }

};
#endif //_TESTMONODOMAINCONDUCTIONVELOCITY_HPP_
