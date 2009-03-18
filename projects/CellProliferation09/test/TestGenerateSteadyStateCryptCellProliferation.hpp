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
#ifndef TESTGENERATESTEEADYSTATECRYPTCELLPROLIFERATION_HPP_
#define TESTGENERATESTEEADYSTATECRYPTCELLPROLIFERATION_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cancer headers
#include "TissueSimulationArchiver.hpp"

#include "CryptSimulation2d.hpp"
#include "StochasticDurationGenerationBasedCellCycleModelCellsGenerator.hpp"
#include "SimpleWntCellCycleModelCellsGenerator.hpp"
#include "StochasticWntCellCycleModelCellsGenerator.hpp"
#include "IngeWntSwatCellCycleModelCellsGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "SloughingCellKiller.hpp"
#include "LogFile.hpp"
#include "SunterSetup.hpp"

#include <cmath>

class TestGenerateSteadyStateCryptCellProliferation : public CxxTest::TestSuite
{
public:

    /*
     * This test generates steady state crypts for use in labelling and monoclonality simulations.
     *
     * The setup is done with:
     * a) the sunter_index (which gives different geometries; 1 and 3 referring to different sites in the mouse colon)
     * b) the cell_cycle_index which takes the following values:
     *    1. For Meineke-style cells with a stochastic element to their cell cycle times.
     *    2. A Simple-Wnt based model which divides when above a certain Wnt threshold.
     *    3. A Stochastic-Wnt model which uses ODEs to decide G1 duration and adds a stochastic element to G2.
     *    4. van Leeuwen Hypothesis 1-Swat cell cycle model.
     *    5. van Leeuwen hypothesis 2-Swat cell cycle model.
     *
     * The archives which are created should be saved and used as an input into other simulations.
     */
    void TestGenerateSteadyStateResultsForCellProliferation09Experiments() throw(Exception)
    {
        unsigned sunter_index = 1u; // the geometry to use (see above)
        unsigned cell_cycle_index = 4u; // the cell cycle model to use (see above)

        // Set up instance of CancerParameters singleton
        CancerParameters* p_params = CancerParameters::Instance();

        // Set output directory
        std::string output_directory = "SteadyStateCryptForLabelling";

        // Set up instance of LogFile singleton
        LogFile::Instance()->Set(1, output_directory, "log.dat");

        double end_of_simulation = 300.0; // hours
        double time_of_each_run = 50.0; // for each run

        SunterSetup sunter_setup(sunter_index); // SunterGeometry and parameters set here
        HoneycombMeshGenerator generator = sunter_setup.SetSunterParametersAndGetMeshGenerator();

        // Create cylindrical mesh
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Get location indices corresponding to real cells in mesh
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);

        // Set up cells
        std::vector<TissueCell> cells;
        switch (cell_cycle_index)
        {
            case 1u:
                {
                    StochasticDurationGenerationBasedCellCycleModelCellsGenerator<2> cell_generator;
                    cell_generator.GenerateForCrypt(cells, *p_mesh, location_indices, true);
                }
                WntConcentration::Instance()->SetType(NONE);
                break;
            case 2u:
                {
                    SimpleWntCellCycleModelCellsGenerator<2> cell_generator;
                    cell_generator.GenerateForCrypt(cells, *p_mesh, location_indices, true);
                }
                WntConcentration::Instance()->SetType(LINEAR);
                break;
            case 3u:
                {
                    StochasticWntCellCycleModelCellsGenerator<2> cell_generator;
                    cell_generator.GenerateForCrypt(cells, *p_mesh, location_indices, true);
                }
                WntConcentration::Instance()->SetType(LINEAR);
                CancerParameters::Instance()->SetTopOfLinearWntConcentration(4.0/18.0);
                break;
            case 4u:
                {
                    IngeWntSwatCellCycleModelCellsGenerator<2> cell_generator(1u);
                    cell_generator.GenerateForCrypt(cells, *p_mesh, location_indices, true);
                }
                WntConcentration::Instance()->SetType(LINEAR);
                CancerParameters::Instance()->SetTopOfLinearWntConcentration(8.0/18.0);
                break;
            case 5u:
                {
                    IngeWntSwatCellCycleModelCellsGenerator<2> cell_generator(2u);
                    cell_generator.GenerateForCrypt(cells, *p_mesh, location_indices, true);
                }
                WntConcentration::Instance()->SetType(LINEAR);
                CancerParameters::Instance()->SetTopOfLinearWntConcentration(11.0/18.0);
                break;
            default:
                EXCEPTION("Please enter a valid index for the cell cycle model");
        }

        // Create crypt facade
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Set up instance of WntConcentration singleton and associate it with crypt
        WntConcentration::Instance()->SetTissue(crypt);

        // Set up force law
        GeneralisedLinearSpringForce<2> meineke_force;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&meineke_force);

        // Create crypt simulation
        CryptSimulation2d simulator(crypt, force_collection);

        // Set where to output simulation results
        simulator.SetOutputDirectory(output_directory);

        // Set simulation to output cell types
        simulator.SetOutputCellTypes(true);

        // Set length of simulation
        simulator.SetEndTime(time_of_each_run);

        // Only save results every tenth time step
        simulator.SetSamplingTimestepMultiple(10);

        // Set up sloughing cell killer and pass in to simulation
        AbstractCellKiller<2>* p_cell_killer = new SloughingCellKiller(&simulator.rGetTissue(), 0.01);
        simulator.AddCellKiller(p_cell_killer);

        // UNUSUAL SET UP HERE /////////////////////////////////////

        p_params->SetDampingConstantNormal(1.0); // normally 1

        // Do not give mutant cells any different movement properties to normal ones
        p_params->SetDampingConstantMutant(p_params->GetDampingConstantNormal());

        p_params->SetSpringStiffness(30.0); // normally 15.0;
        // 0.3/30 = 0.01 (i.e. Meineke's values)

        simulator.UseJiggledBottomCells();

        // END OF UNUSUAL SET UP //////////////////////////////////

        // Run simulation
        simulator.Solve();

        // Save simulation
        TissueSimulationArchiver<2, CryptSimulation2d>::Save(&simulator);

        for (double t=time_of_each_run; t<end_of_simulation; t+=time_of_each_run)
        {
            // Load simulation
            CryptSimulation2d* p_simulator = TissueSimulationArchiver<2, CryptSimulation2d>::Load(output_directory, t);

            // Reset end time and run simulation
            p_simulator->SetEndTime(t+time_of_each_run);
            p_simulator->Solve();

            // Save and tidy up
            TissueSimulationArchiver<2, CryptSimulation2d>::Save(p_simulator);
            delete p_simulator;
        }

        // Close log file and tidy up
        LogFile::Close();
        delete p_cell_killer;
        delete p_params;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        WntConcentration::Destroy();
    }
};

#endif /*TESTGENERATESTEEADYSTATECRYPTCELLPROLIFERATION_HPP_*/
