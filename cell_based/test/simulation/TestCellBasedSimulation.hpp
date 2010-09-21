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
#ifndef TESTCELLBASEDSIMULATION_HPP_
#define TESTCELLBASEDSIMULATION_HPP_

#include <cxxtest/TestSuite.h>

#include <cstdio>
#include <ctime>
#include <cmath>

#include "CellBasedSimulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CryptCellsGenerator.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "LinearSpringWithVariableSpringConstantsForce.hpp"
#include "RandomCellKiller.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NumericFileComparison.hpp"
#include "CellBasedEventHandler.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StochasticWntCellCycleModel.hpp"

// Simple subclass of CellBasedSimulation which just overloads StoppingEventHasOccurred
// for testing the stopping event functionality..
class CellBasedSimulationWithMyStoppingEvent : public CellBasedSimulation<2>
{
private:
    // define a stopping event with says stop if t>3.14
    bool StoppingEventHasOccurred()
    {
        return  (SimulationTime::Instance()->GetTime() > 3.1415);
    }

public:
    CellBasedSimulationWithMyStoppingEvent(AbstractCellPopulation<2>& rCellPopulation,
                                        std::vector<AbstractForce<2>* > forceCollection)
      : CellBasedSimulation<2>(rCellPopulation, forceCollection)
   {
   }
};


/**
 *  Note: Most tests of CellBasedSimulation are in TestCryptSimulation2d
 */
class TestCellBasedSimulation : public AbstractCellBasedTestSuite
{
private:

    double mLastStartTime;
    void setUp()
    {
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

    void TestOutputStatistics() throw(Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        // Set up mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<StochasticWntCellCycleModel> cell_generator;
        cell_generator.Generate(cells, p_mesh, location_indices, true);

        // Set up cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetOutputCellPopulationVolumes(true);
        cell_population.SetOutputCellVariables(true);
        cell_population.SetOutputCellCyclePhases(true);
        cell_population.SetOutputCellAges(true);

        // Set up Wnt Gradient
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(cell_population);

        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.UseCutoffPoint(1.5);
        std::vector<AbstractForce<2>* > force_collection;
        force_collection.push_back(&linear_force);

        // Set up cell-based simulation
        CellBasedSimulation<2> simulator(cell_population, force_collection);
        simulator.SetOutputDirectory("CellBasedSimulationWritingProteins");
        simulator.SetEndTime(0.5);

        TS_ASSERT_DELTA(simulator.GetDt(), 1.0/120.0, 1e-12);

        // Run cell-based simulation
        TS_ASSERT_EQUALS(simulator.GetOutputDirectory(), "CellBasedSimulationWritingProteins");
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        OutputFileHandler handler("CellBasedSimulationWritingProteins", false);

        std::string cell_variables_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/cellvariables.dat";
        NumericFileComparison comp_cell_variables(cell_variables_file, "cell_based/test/data/CellBasedSimulationWritingProteins/cellvariables.dat");
        TS_ASSERT(comp_cell_variables.CompareFiles(1e-2));

        std::string cell_cycle_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/cellcyclephases.dat";
        NumericFileComparison comp_cell_cycle(cell_cycle_file, "cell_based/test/data/CellBasedSimulationWritingProteins/cellcyclephases.dat");
        TS_ASSERT(comp_cell_cycle.CompareFiles(1e-2));

        std::string cell_ages_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/cellages.dat";
        NumericFileComparison comp_cell_ages(cell_ages_file, "cell_based/test/data/CellBasedSimulationWritingProteins/cellages.dat");
        TS_ASSERT(comp_cell_ages.CompareFiles(1e-2));

        // Tidy up
        WntConcentration<2>::Destroy();
    }

    void TestOutputNodeVelocities() throw(Exception)
    {
        // Create a simple mesh
        HoneycombMeshGenerator generator(5, 5, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Create a differentiated cell for each node
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(DIFFERENTIATED);

            double birth_time = -RandomNumberGenerator::Instance()->ranf()*
                                (CellBasedConfig::Instance()->GetStemCellG1Duration()
                                    + CellBasedConfig::Instance()->GetSG2MDuration() );

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create a cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create a force law
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.UseCutoffPoint(1.5);
        std::vector<AbstractForce<2>* > force_collection;
        force_collection.push_back(&linear_force);

        // Set up simulation
        CellBasedSimulation<2> simulator(cell_population, force_collection);
        simulator.SetOutputDirectory("TestOutputNodeVelocities");
        simulator.SetEndTime(0.5);

        // Record node velocities
        TS_ASSERT_EQUALS(simulator.GetOutputNodeVelocities(), false);
        simulator.SetOutputNodeVelocities(true);

        // Run simulation
        simulator.Solve();

        // Check node velocities file
        // The velocities should all be zero(ish), as the cell population is in mechanical equilibrium
        OutputFileHandler handler("TestOutputNodeVelocities", false);

        std::string node_velocities_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/nodevelocities.dat";
        NumericFileComparison node_velocities(node_velocities_file, "cell_based/test/data/TestOutputNodeVelocities/nodevelocities.dat");
        TS_ASSERT(node_velocities.CompareFiles(1e-2));
    }

    void TestOutputNodeVelocitiesWithGhostNodes() throw(Exception)
    {
        // Create a simple mesh with a surrounding layer of ghost nodes
        HoneycombMeshGenerator generator(3, 3, 1, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create a differentiated cell for each non-ghost node
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<location_indices.size(); i++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(DIFFERENTIATED);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetBirthTime(-1.0);

            cells.push_back(p_cell);
        }

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        // Create a force law
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.UseCutoffPoint(1.5);
        std::vector<AbstractForce<2>* > force_collection;
        force_collection.push_back(&linear_force);

        // Set up simulation
        CellBasedSimulation<2> simulator(cell_population, force_collection);
        simulator.SetOutputDirectory("TestOutputNodeVelocitiesWithGhostNodes");
        simulator.SetEndTime(0.5);

        // Record node velocities
        simulator.SetOutputNodeVelocities(true);

        // Run simulation
        simulator.Solve();
    }

    /**
     * Test a cell-based simulation with a cell killer.
     *
     * In this test, we solve a cell-based simulation without ghost nodes and
     * check that the numbers of nodes and cells match at the end of the
     * simulation.
     */
    void TestCellBasedSimulationWithCellDeath() throw (Exception)
    {
        // Create a simple mesh
        int num_cells_depth = 5;
        int num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells, one for each node. Give each cell a random birth time.
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(STEM);

            double birth_time = -RandomNumberGenerator::Instance()->ranf()*
                                (CellBasedConfig::Instance()->GetStemCellG1Duration()
                                    + CellBasedConfig::Instance()->GetSG2MDuration() );
            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create a cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create a force law
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.UseCutoffPoint(1.5);
        std::vector<AbstractForce<2>* > force_collection;
        force_collection.push_back(&linear_force);

        // Set up cell-based simulation
        CellBasedSimulation<2> simulator(cell_population, force_collection);
        simulator.SetOutputDirectory("TestCellBasedSimulationWithCellDeath");
        simulator.SetEndTime(0.5);

        // Add cell killer
        RandomCellKiller<2> random_cell_killer(&cell_population, 0.997877574);
        simulator.AddCellKiller(&random_cell_killer);

        // For coverage of an exception.
        simulator.SetUpdateCellPopulationRule(false);
        TS_ASSERT_THROWS_THIS(simulator.Solve(),"CellPopulation has had births or deaths but mUpdateCellPopulation is set to false, please set it to true.");
        CellBasedEventHandler::Reset(); // Otherwise logging has been started but not stopped due to exception above.

        simulator.SetUpdateCellPopulationRule(true);
        simulator.Solve();

        // Check that the number of nodes is equal to the number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumNodes(), simulator.rGetCellPopulation().GetNumRealCells());

        // For coverage of these 'Get' functions
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 17u);
    }

    /**
	 * Test a cell-based simulation with multiple forces.
	 */
	void TestCellBasedSimulationWithMultipleForces() throw (Exception)
	{
		// Create a simple mesh
		int num_cells_depth = 5;
		int num_cells_width = 5;
		HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
		MutableMesh<2,2>* p_mesh = generator.GetMesh();

		// Set up cells, one for each node. Give each cell a random birth time.
		std::vector<CellPtr> cells;
		boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
		for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
		{
			FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
			p_model->SetCellProliferativeType(STEM);

			double birth_time = -RandomNumberGenerator::Instance()->ranf()*
								(CellBasedConfig::Instance()->GetStemCellG1Duration()
									+ CellBasedConfig::Instance()->GetSG2MDuration() );
			CellPtr p_cell(new Cell(p_state, p_model));
			p_cell->SetBirthTime(birth_time);
			cells.push_back(p_cell);
		}

		// Create a cell population
		MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

		// Create a force law system
		std::vector<AbstractForce<2>* > force_collection;
		GeneralisedLinearSpringForce<2> linear_force;
		force_collection.push_back(&linear_force);
		LinearSpringWithVariableSpringConstantsForce<2> variable_force;
		force_collection.push_back(&variable_force);

		// Set up cell-based simulation
		CellBasedSimulation<2> simulator(cell_population, force_collection);
		simulator.SetOutputDirectory("TestCellBasedSimulationWithMultipleForces");
		simulator.SetEndTime(0.5);

		simulator.Solve();

		// Check that the number of nodes is equal to the number of cells
		TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumNodes(), simulator.rGetCellPopulation().GetNumRealCells());

		// For coverage of these 'Get' functions
		TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
		TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
	}


    void TestCellBasedSimulationWithStoppingEvent() throw (Exception)
    {
        HoneycombMeshGenerator generator(2, 2, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells, one for each node. Give each cell a random birth time.
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        // Create a cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create a force law
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.UseCutoffPoint(1.5);
        std::vector<AbstractForce<2>* > force_collection;
        force_collection.push_back(&linear_force);

        // Set up cell-based simulation WITH the stopping event
        CellBasedSimulationWithMyStoppingEvent simulator(cell_population, force_collection);
        simulator.SetOutputDirectory("TestCellPopulationSimWithStoppingEvent");

        // Set the end time to 10.0 - the stopping event is, however, t>3.1415.
        simulator.SetEndTime(10.0);

        // Run cell-based simulation
        simulator.Solve();

        double time = SimulationTime::Instance()->GetTime();
        TS_ASSERT_DELTA(time, 3.1415, 1e-1); // big tol, doesn't matter, just want t~3.14 and t!=10
        // t should be strictly greater than the 3.1415
        TS_ASSERT_LESS_THAN(3.1415, time);
    }


    void TestApoptosisSpringLengths() throw (Exception)
    {
        unsigned num_cells_depth = 2;
        unsigned num_cells_width = 2;
        double crypt_length = num_cells_depth-0.0;
        double crypt_width = num_cells_width-0.0;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        CellBasedConfig* p_params = CellBasedConfig::Instance();
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);

        // Set up cells
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<location_indices.size(); i++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(TRANSIT);

            CellPtr p_cell(new Cell(p_state, p_model));

            double birth_time = -RandomNumberGenerator::Instance()->ranf()*(p_params->GetTransitCellG1Duration()
                                               +p_params->GetSG2MDuration());

            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        // Create force law
        GeneralisedLinearSpringForce<2> linear_force;
        std::vector<AbstractForce<2>* > force_collection;
        force_collection.push_back(&linear_force);

        // Create crypt simulation from cell population and force law
        CellBasedSimulation<2> simulator(cell_population, force_collection);
        simulator.SetOutputDirectory("2dSpheroidApoptosis");
        simulator.SetEndTime(1.0);

        cell_population.GetCellUsingLocationIndex(14)->SetApoptosisTime(2.0);
        cell_population.GetCellUsingLocationIndex(15)->SetApoptosisTime(2.0);
        cell_population.GetCellUsingLocationIndex(14)->StartApoptosis();
        cell_population.GetCellUsingLocationIndex(15)->StartApoptosis();
        simulator.SetNoBirth(true);

        // Run cell-based simulation
        simulator.Solve();

        /*
         * We track the locations of two dying cells (a and b) and two
         * live cells adjacent to them (c and d)
         *
         * All cells begin distance 1 apart.
         *
         * a and b move together to leave a gap of 0.
         * a and c (and b and d) move to a distance of 0.5 apart.
         */

        c_vector<double, 2> a_location = cell_population.rGetMesh().GetNode(14)->rGetLocation();
        c_vector<double, 2> b_location = cell_population.rGetMesh().GetNode(15)->rGetLocation();
        c_vector<double, 2> c_location = cell_population.rGetMesh().GetNode(20)->rGetLocation();
        c_vector<double, 2> d_location = cell_population.rGetMesh().GetNode(21)->rGetLocation();

        double a_b_separation = sqrt((a_location[0]-b_location[0])*(a_location[0]-b_location[0]) +
                                (a_location[1]-b_location[1])*(a_location[1]-b_location[1]));
        double a_c_separation = sqrt((a_location[0]-c_location[0])*(a_location[0]-c_location[0]) +
                                (a_location[1]-c_location[1])*(a_location[1]-c_location[1]));
        double c_d_separation = sqrt((d_location[0]-c_location[0])*(d_location[0]-c_location[0]) +
                                (d_location[1]-c_location[1])*(d_location[1]-c_location[1]));

        TS_ASSERT_DELTA(a_b_separation, 0.5, 1e-1);
        TS_ASSERT_DELTA(a_c_separation, 0.75, 1e-1);
        TS_ASSERT_DELTA(c_d_separation, 1.0, 1e-1);

        // Reset end time and run cell-based simulation
        simulator.SetEndTime(1.99);
        simulator.Solve();

        a_location = cell_population.rGetMesh().GetNode(14)->rGetLocation();
        b_location = cell_population.rGetMesh().GetNode(15)->rGetLocation();
        c_location = cell_population.rGetMesh().GetNode(20)->rGetLocation();
        d_location = cell_population.rGetMesh().GetNode(21)->rGetLocation();

        a_b_separation = sqrt((a_location[0]-b_location[0])*(a_location[0]-b_location[0]) +
                         (a_location[1]-b_location[1])*(a_location[1]-b_location[1]));
        a_c_separation = sqrt((a_location[0]-c_location[0])*(a_location[0]-c_location[0]) +
                         (a_location[1]-c_location[1])*(a_location[1]-c_location[1]));
        c_d_separation = sqrt((d_location[0]-c_location[0])*(d_location[0]-c_location[0]) +
                         (d_location[1]-c_location[1])*(d_location[1]-c_location[1]));

        TS_ASSERT_DELTA(a_b_separation, 0.01, 1e-1);
        TS_ASSERT_DELTA(a_c_separation, 0.5, 1e-1);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 4u);

        // Reset end time and run cell-based simulation
        simulator.SetEndTime(2.01);
        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 2u);
    }


    void TestCellBasedSimulationParameterOutputMethods() throw (Exception)
    {
        // Create a simple mesh
        int num_cells_depth = 5;
        int num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        // Create a cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create a force law Collection
        std::vector<AbstractForce<2>* > force_collection;

        // Set up cell-based simulation
        CellBasedSimulation<2> simulator(cell_population, force_collection);

        ///\todo uncomment see #1453
        //TS_ASSERT_EQUALS(simulator.GetIdentifier(), "CellBasedSimulation<2>");

		std::string output_directory = "TestCellBasedSimulationOutputParameters";
		OutputFileHandler output_file_handler(output_directory, false);
		out_stream parameter_file = output_file_handler.OpenOutputFile("cell_based_sim_results.parameters");
		simulator.OutputSimulationParameters(parameter_file);
		parameter_file->close();

		std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();
		TS_ASSERT_EQUALS(system(("diff " + results_dir + "cell_based_sim_results.parameters  cell_based/test/data/TestCellBasedSimulationOutputParameters/cell_based_sim_results.parameters").c_str()), 0);

		simulator.SetOutputDirectory("TestCellBasedSimulationOutputParameters");
		simulator.OutputSimulationSetup();
		///\todo #1453 This is to do with the pre Boost 1.37 problem --- 
        //TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.parameters  cell_based/test/data/TestCellBasedSimulationOutputParameters/results.parameters").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff --ignore-matching-lines=\"CellPopulation\" " + results_dir + "results.parameters  cell_based/test/data/TestCellBasedSimulationOutputParameters/results.parameters").c_str()), 0);

		//Check that the files which we don't want to compare actually exist
		std::ifstream machine_file;
		std::string command = results_dir+"/system.info.0";
		machine_file.open(command.c_str());
		TS_ASSERT(machine_file.is_open());
		machine_file.close();

		std::ifstream info_file;
		command = results_dir+"/build.info";
		info_file.open(command.c_str());
		TS_ASSERT(info_file.is_open());
		info_file.close();
    }
};
#endif /*TESTCELLBASEDSIMULATION_HPP_*/
