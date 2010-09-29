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
#ifndef TESTCRYPTSIMULATION2DNIGHTLY_HPP_
#define TESTCRYPTSIMULATION2DNIGHTLY_HPP_

#include <cxxtest/TestSuite.h>

#include "CryptSimulation2d.hpp"
#include "CryptCellsGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "WntConcentration.hpp"
#include "RandomCellKiller.hpp"
#include "SloughingCellKiller.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "CellBasedEventHandler.hpp"
#include "NumericFileComparison.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"


#include "WntCellCycleModel.hpp"
#include "SimpleWntCellCycleModel.hpp"
class TestCryptSimulation2dNightly : public AbstractCellBasedTestSuite
{
private:

    double mLastStartTime;
    void setUp()
    {
        CellBasedEventHandler::Disable(); // these tests fail with event-handling on
        mLastStartTime = std::clock();
        AbstractCellBasedTestSuite::setUp();
    }

    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCellBasedTestSuite::tearDown();
    }

public:
    /**
     * Sloughing with a sloughing cell killer and not turning
     * into ghost nodes on a non-periodic mesh
     */
    void TestSloughingCellKillerOnNonPeriodicCrypt() throw (Exception)
    {
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;

        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create force law
        GeneralisedLinearSpringForce<2> linear_force;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Create crypt simulation from cell population and force law
        CryptSimulation2d simulator(crypt, force_collection);
        simulator.SetOutputDirectory("Crypt2DSloughingDeathNonPeriodic");
        simulator.SetEndTime(4.0);

        // Create cell killer and pass in to crypt simulation
        SloughingCellKiller<2> sloughing_cell_killer(&crypt,
                                                     CellBasedConfig::Instance()->GetCryptLength()
                                                     true,
                                                     CellBasedConfig::Instance()->GetCryptWidth());

        simulator.AddCellKiller(&sloughing_cell_killer);

        // Run simulation
        simulator.Solve();
    }

    void TestSloughingDeathWithPeriodicMesh() throw (Exception)
    {
        unsigned cells_across = 7;
        unsigned cells_up = 12;
        double crypt_width = 6.0;
        unsigned thickness_of_ghost_layer = 4;

        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer,true,crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create force law
        GeneralisedLinearSpringForce<2> linear_force;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Create crypt simulation from cell population and force law
        CryptSimulation2d simulator(crypt, force_collection);
        simulator.SetOutputDirectory("Crypt2DSloughingDeathPeriodic");
        simulator.SetEndTime(4.0);

        // Create cell killer and pass in to crypt simulation
        SloughingCellKiller<2> cell_killer(&crypt, CellBasedConfig::Instance()->GetCryptLength());
        simulator.AddCellKiller(&cell_killer);

        // Run simulation
        simulator.Solve();

        std::vector<bool> ghost_node_indices_after = (static_cast<MeshBasedCellPopulationWithGhostNodes<2>* >(&(simulator.rGetCellPopulation())))->rGetGhostNodes();
        unsigned num_ghosts = 0;
        for (unsigned i=0; i<ghost_node_indices_after.size(); i++)
        {
            if (ghost_node_indices_after[i])
            {
                num_ghosts++;
            }
        }

        // Check no new ghost nodes have been created
        TS_ASSERT_EQUALS(num_ghosts, p_mesh->GetNumNodes() - crypt.GetNumRealCells());

        // There should be this number of cells left after this amount of time
        // (we have lost two rows of 7 but had a bit of birth too)
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 85u);
    }



    void TestMonolayerWithCutoffPointAndNoGhosts() throw (Exception)
    {
        unsigned num_cells_depth = 11;
        unsigned num_cells_width = 6;
        double crypt_length = num_cells_depth-1.0;
        double crypt_width = num_cells_width-1.0;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        CellBasedConfig::Instance()->SetCryptLength(crypt_length);
        CellBasedConfig::Instance()->SetCryptWidth(crypt_width);

        // Set up cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true, -1.0);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create force law
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.SetCutOffLength(sqrt(2)); // root2 is a sensible choice
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Create crypt simulation from cell population and force law
        CryptSimulation2d simulator(crypt, force_collection);
        simulator.SetOutputDirectory("MonolayerCutoffPointNoGhosts");
        simulator.SetEndTime(12.0);

        // Run simulation
        simulator.Solve();
    }

    /*
	* This tests that the results files are correct (added because of #1130).
	*/
	void TestResultsFileForLongerCryptSimulation() throw(Exception)
	{
		// Set output directory
		std::string output_directory = "TestResultsFileForLongerCryptSimulation";

		// Create cylindrical mesh
		HoneycombMeshGenerator generator(16, 19, 0, true);
		Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

		// Get location indices corresponding to real cells in mesh
		std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

		// Set up each cell with a simple Wnt-based cell cycle model
		std::vector<CellPtr> cells;
		CryptCellsGenerator<SimpleWntCellCycleModel> cell_generator;
		cell_generator.Generate(cells, p_mesh, location_indices, true);

		// Set some model parameters for the cell cycle model
		for(unsigned index=0; index < cells.size(); index++)
		{
		   cells[index]->GetCellCycleModel()->SetSDuration(7.4);
		   cells[index]->GetCellCycleModel()->SetG2Duration(1.4);
		   cells[index]->GetCellCycleModel()->SetMDuration(0.72);
		   cells[index]->GetCellCycleModel()->SetTransitCellG1Duration(9.4);
		   cells[index]->GetCellCycleModel()->SetStemCellG1Duration(9.4);
		}

		// Create cell population
		MeshBasedCellPopulation<2> crypt(*p_mesh, cells);

		// Set cell population to output cell types
		crypt.SetOutputCellProliferativeTypes(true);

		// Set up instance of WntConcentration singleton and associate it with crypt
		WntConcentration<2>::Instance()->SetType(LINEAR);
		WntConcentration<2>::Instance()->SetCellPopulation(crypt);

		// Set up force law
		GeneralisedLinearSpringForce<2> meineke_force;

		// Unusual set-up here (corresponds to the Meineke crypt model parameters)
		meineke_force.SetMeinekeSpringStiffness(30.0);
		// Sets the MeinekeSpringGrowthDuration to be the default MPhase Duration
		meineke_force.SetMeinekeSpringGrowthDuration(cells[0]->GetCellCycleModel()->GetMDuration());

		// Pass force law into collection
		std::vector<AbstractForce<2>*> force_collection;
		force_collection.push_back(&meineke_force);

		// Create crypt simulation
		CryptSimulation2d simulator(crypt, force_collection);

		// Set where to output simulation results
		simulator.SetOutputDirectory(output_directory);

		// Set length of simulation
		simulator.SetEndTime(20.0);

		// Only save results every tenth time step
		simulator.SetSamplingTimestepMultiple(10);

		// Set up sloughing cell killer and pass in to simulation
		AbstractCellKiller<2>* p_cell_killer = new SloughingCellKiller<2>(&simulator.rGetCellPopulation(), CellBasedConfig::Instance()->GetCryptLength());
		simulator.AddCellKiller(p_cell_killer);

		// Unusual set-up here (corresponds to the Meineke crypt model parameters)
		simulator.UseJiggledBottomCells();

		// Run simulation
		simulator.Solve();

		// Test that results files are correct
		OutputFileHandler handler(output_directory, false);
		std::string results_dir = handler.GetOutputDirectoryFullPath() + "results_from_time_0";

		NumericFileComparison comp_ele(results_dir + "/results.vizelements", "cell_based/test/data/TestResultsFileForLongerCryptSimulation/results.vizelements");
		TS_ASSERT(comp_ele.CompareFiles());
		TS_ASSERT_EQUALS(system(("diff " + results_dir + "/results.vizelements cell_based/test/data/TestResultsFileForLongerCryptSimulation/results.vizelements").c_str()), 0);

		NumericFileComparison comp_nodes(results_dir + "/results.viznodes", "cell_based/test/data/TestResultsFileForLongerCryptSimulation/results.viznodes");
		TS_ASSERT(comp_nodes.CompareFiles(1e-15));

		NumericFileComparison comp_celltypes(results_dir + "/results.vizcelltypes", "cell_based/test/data/TestResultsFileForLongerCryptSimulation/results.vizcelltypes");
		TS_ASSERT(comp_celltypes.CompareFiles(1e-15));

		TS_ASSERT_EQUALS(system(("diff " + results_dir + "/results.vizsetup cell_based/test/data/TestResultsFileForLongerCryptSimulation/results.vizsetup").c_str()), 0);

		// Tidy up
		delete p_cell_killer;
		WntConcentration<2>::Destroy();
	}
};


#endif /*TESTCRYPTSIMULATION2DNIGHTLY_HPP_*/
