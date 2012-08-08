/*

Copyright (c) 2005-2012, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/
#ifndef TESTGENERATESTEADYSTATECRYPTFORINVASION_HPP_
#define TESTGENERATESTEADYSTATECRYPTFORINVASION_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "CryptSimulation2dWithCryptInvasionStoppingEvent.hpp"
#include "CryptCellsGenerator.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "SloughingCellKiller.hpp"
#include "SunterSetup.hpp"
#include "SmartPointers.hpp"

#include <cmath>

class TestGenerateSteadyStateCryptForInvasion : public CxxTest::TestSuite
{
public:

    /**
     * This test generates steady state crypts for use in crypt invasion simulations.
     *
     * The setup is done with the sunter_index (which gives different geometries;
     * 1 and 3 referring to different sites in the mouse colon)
     *
     * The archives which are created should be saved and used as an input into other simulations.
     */
    void TestGenerateSteadyStateResultsForCryptInvasionExperiments() throw(Exception)
    {
        std::cout << "\nGenerating steady state archives...\n" << std::flush;

        unsigned sunter_index = 3u; // the geometry to use (see above)

        // Set output directory
        std::string output_directory = "SteadyStateCryptForInvasion";

        double end_of_simulation = 300.0; // hours
        double time_of_each_run = 50.0; // for each run

        SunterSetup sunter_setup(sunter_index); // SunterGeometry and parameters set here
        CylindricalHoneycombMeshGenerator generator = sunter_setup.SetSunterParametersAndGetMeshGenerator();

        // Store the height of the crypt, for use by the SloughingCellKiller
        double crypt_height = generator.GetDomainDepth();

        // Create cylindrical mesh
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Get location indices corresponding to real cells in mesh
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up instance of SimulationTime singleton
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);

        // Set up each cell with a simple Wnt-based cell cycle model
        std::vector<boost::shared_ptr<Cell> > cells;
        CryptCellsGenerator<SimpleWntCellCycleModel> cell_generator;
        cell_generator.Generate(cells, p_mesh, location_indices, true);
        sunter_setup.SetUpCellCycleModelParameters(cells);

        // Create crypt facade
        MeshBasedCellPopulation<2> crypt(*p_mesh, cells);

        // Set up instance of WntConcentration singleton and associate it with crypt
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_height);

        // Create crypt simulation
        CryptSimulation2dWithCryptInvasionStoppingEvent simulator(crypt, crypt_height);

        // Set where to output simulation results
        simulator.SetOutputDirectory(output_directory);

        // Set up force law
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_meineke_force);
        p_meineke_force->SetMeinekeSpringStiffness(30.0);
        simulator.AddForce(p_meineke_force);

        // Set simulation to output cell types
        crypt.SetOutputCellProliferativeTypes(true);
        crypt.SetOutputCellMutationStates(true);

        // Set length of simulation
        simulator.SetEndTime(time_of_each_run);

        // Only write results to file every hour
        simulator.SetSamplingTimestepMultiple(120);

        // Set up sloughing cell killer and pass in to simulation
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_cell_killer, (&simulator.rGetCellPopulation(), crypt_height));
        simulator.AddCellKiller(p_cell_killer);

        // UNUSUAL SET UP HERE /////////////////////////////////////

        crypt.SetDampingConstantNormal(1.0); // normally 1

        // Do not give mutant cells any different movement properties to normal ones
        crypt.SetDampingConstantMutant(crypt.GetDampingConstantNormal());

        // A small random upward force is applied to cells at the base to
        // prevent an unstable equilibrium with many cells with compressed springs
        // building up on y=0.
        
        simulator.UseJiggledBottomCells();
		
		// Don't check for stopping events whilst getting to 
		// a steady state.
		simulator.SetCheckForStoppingEvent(false);
		
        // END OF UNUSUAL SET UP //////////////////////////////////

        // Run simulation
        simulator.Solve();

        // Save simulation
        CellBasedSimulationArchiver<2, CryptSimulation2dWithCryptInvasionStoppingEvent>::Save(&simulator);

        unsigned count = 0;
        unsigned num_cells;
        unsigned num_wild_type_cells;
        for (double t=time_of_each_run; t<end_of_simulation; t+=time_of_each_run)
        {
            // Load simulation
            CryptSimulation2dWithCryptInvasionStoppingEvent* p_simulator =
                    CellBasedSimulationArchiver<2, CryptSimulation2dWithCryptInvasionStoppingEvent>::Load(output_directory, t);

            // Reset end time and run simulation
            p_simulator->SetEndTime(t+time_of_each_run);
            p_simulator->Solve();

            std::vector<boost::shared_ptr<Cell> > bottom_row = p_simulator->GetBottomRowOfCells();
            TS_ASSERT_LESS_THAN(bottom_row.size(), 22u);
            count++;

            bool mixed_bottom_row = p_simulator->VectorContainsNormalAndMutantCells(bottom_row);
            TS_ASSERT_EQUALS(mixed_bottom_row, false);
            TS_ASSERT_THROWS_THIS(p_simulator->GetAverageVerticalForces(bottom_row),
                                 "Need a mixed population of wild-type and mutant cells for this method.");
            TS_ASSERT_EQUALS(p_simulator->GetTimeOverWhichForcesAveraged(), 0u);

            // Save and tidy up
            CellBasedSimulationArchiver<2, CryptSimulation2dWithCryptInvasionStoppingEvent>::Save(p_simulator);

            num_cells = p_simulator->rGetCellPopulation().GetNumRealCells();
            num_wild_type_cells = p_simulator->rGetCellPopulation().GetCellPropertyRegistry()->Get<WildTypeCellMutationState>()->GetCellCount();
            delete p_simulator;
        }

        // Check mutation states are reported properly.
        // NB. Labelled cells are no longer "mutation states".
        crypt.WriteResultsToFiles(); // Update to last time point
        std::vector<unsigned> mutation_states_after = crypt.GetCellMutationStateCount();
        TS_ASSERT_EQUALS(mutation_states_after.size(), 4u);
        TS_ASSERT_EQUALS(mutation_states_after[0], crypt.GetNumRealCells());
        TS_ASSERT_EQUALS(mutation_states_after[1], 0u);
        TS_ASSERT_EQUALS(mutation_states_after[2], 0u);
        TS_ASSERT_EQUALS(mutation_states_after[3], 0u);
        TS_ASSERT_EQUALS(num_cells,num_wild_type_cells);

        // Close log file and tidy up
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        WntConcentration<2>::Destroy();
    }
};

#endif /*TESTGENERATESTEADYSTATECRYPTFORINVASION_HPP_*/
