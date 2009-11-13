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
#ifndef TESTGENERATESTEADYSTATECRYPT_HPP_
#define TESTGENERATESTEADYSTATECRYPT_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cancer headers
#include "TissueSimulationArchiver.hpp"

#include "CryptSimulation2d.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "StochasticWntCellCycleModelCellsGenerator.hpp"
#include "SloughingCellKiller.hpp"


class TestGenerateSteadyStateCrypt : public CxxTest::TestSuite
{
public:

    /*
     * This test can be used to generate steady state crypts for use
     * as the starting points of other simulations.
     *
     * You need to specify :
     * the kind of cell cycle model to use on line 64,
     * WntConcentration on line 69,
     * change any cancer parameters around line 90,
     * and give the simulator options around line 95.
     */
    void TestGenerateSteadyStateCryptArchives() throw (Exception)
    {
        TissueConfig* p_params = TissueConfig::Instance();
        std::string output_directory = "SteadyStateCrypt";

        double end_of_simulation = 150.0; // hours
        double time_of_each_run = 10.0; // for each run - the more saves and loads the better for testing this

        // Create mesh
        unsigned cells_across = 13;
        unsigned cells_up = 25;
        double crypt_width = 12.1;
        unsigned thickness_of_ghost_layer = 3;

        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, true, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);

        // Set up cells
        std::vector<TissueCell> cells;
        StochasticWntCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, location_indices, true);

        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        WntConcentration<2>::Instance()->SetType(LINEAR);
        TissueConfig::Instance()->SetTopOfLinearWntConcentration(1.0/3.0);
        WntConcentration<2>::Instance()->SetTissue(crypt);

        GeneralisedLinearSpringForce<2> linear_force;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        CryptSimulation2d simulator(crypt, force_collection);
        simulator.SetOutputDirectory(output_directory);

        // Set simulation to output cell types
        TissueConfig::Instance()->SetOutputCellMutationStates(true);

        // Set length of simulation here
        simulator.SetEndTime(time_of_each_run);

        SloughingCellKiller<2> cell_killer(&simulator.rGetTissue(),0.01);
        simulator.AddCellKiller(&cell_killer);

        // UNUSUAL SET UP HERE /////////////////////////////////////

        p_params->SetDampingConstantNormal(1.0);    // normally 1

        // Do not give mutant cells any different movement properties to normal ones
        p_params->SetDampingConstantMutant(p_params->GetDampingConstantNormal());

        p_params->SetSpringStiffness(30.0); //normally 15.0;
        // 0.3/30 = 0.01 (i.e. Meineke's values)

        simulator.UseJiggledBottomCells();

        // END OF UNUSUAL SET UP! //////////////////////////////////

        simulator.Solve();
        TissueSimulationArchiver<2, CryptSimulation2d>::Save(&simulator);

        for (double t=time_of_each_run; t<end_of_simulation+0.5; t += time_of_each_run)
        {
            CryptSimulation2d* p_simulator = TissueSimulationArchiver<2, CryptSimulation2d>::Load("SteadyStateCrypt",t);
            p_simulator->SetEndTime(t+time_of_each_run);
            p_simulator->Solve();
            TissueSimulationArchiver<2, CryptSimulation2d>::Save(p_simulator);
            delete p_simulator;
        }

        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        WntConcentration<2>::Destroy();
    }

};

#endif /*TESTGENERATESTEADYSTATECRYPT_HPP_*/
