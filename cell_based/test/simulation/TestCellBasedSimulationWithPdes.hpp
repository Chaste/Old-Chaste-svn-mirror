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
#ifndef TESTCELLBASEDSIMULATIONWITHPDES_HPP_
#define TESTCELLBASEDSIMULATIONWITHPDES_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "CellBasedSimulationWithPdes.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "OxygenBasedCellKiller.hpp"
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "SimpleUniformSourcePde.hpp"
#include "CellwiseSourcePde.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "ReplicatableVector.hpp"
#include "NumericFileComparison.hpp"
#include "WildTypeCellMutationState.hpp"

class SimplePdeForTesting : public AbstractLinearEllipticPde<2,2>
{
public:
    double ComputeConstantInUSourceTerm(const ChastePoint<2>& x)
    {
        return -1.0;
    }

    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<2>& x, Element<2,2>*)
    {
        return 0.0;
    }

    c_matrix<double,2,2> ComputeDiffusionTerm(const ChastePoint<2>& )
    {
        return identity_matrix<double>(2);
    }
};


class TestCellBasedSimulationWithPdes : public AbstractCellBasedTestSuite
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

    /*
     * A two-part test for the PostSolve() method.
     *
     * Firstly, test the PDE solver using the problem del squared C = 0.1
     * on the unit disc, with boundary condition C=1 on r=1, which has
     * analytic solution C = 1-0.25*(1-r^2).
     *
     * Secondly, test that cells' hypoxic durations are correctly updated when a
     * nutrient distribution is prescribed.
     */
    void TestPostSolve() throw(Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        // Set up mesh
        MutableMesh<2,2> mesh;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(STEM);
            p_model->SetHypoxicConcentration(0.9);
            p_model->SetQuiescentConcentration(0.9);

            //Use non default G1Durations
			p_model->SetStemCellG1Duration(8.0);
			p_model->SetTransitCellG1Duration(8.0);


            CellPtr p_cell(new Cell(p_state, p_model));

            double birth_time = -RandomNumberGenerator::Instance()->ranf()*18.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Set up cell population
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Set up CellwiseData and associate it with the cell population
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);

        // Since values are first passed in to CellwiseData before it is updated in PostSolve(),
        // we need to pass it some initial conditions to avoid memory errors
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, mesh.GetNode(i)->GetIndex());
        }

        // Set up PDE
        SimplePdeForTesting pde;
        double boundary_value = 1.0;
        bool is_neumann_bc = false;
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, boundary_value, is_neumann_bc);
        std::vector<PdeAndBoundaryConditions<2>*> pde_and_bc_collection;
        pde_and_bc_collection.push_back(&pde_and_bc);

        // Set up force law
        GeneralisedLinearSpringForce<2> linear_force;
        // Use an extremely small cutoff so that no cells interact
        // - this is to ensure that in the Solve method, the cells don't move
        // (we need to call Solve to set up the .vizpdesolution file)
        linear_force.SetCutOffLength(0.0001);
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Set up cell-based simulation
        CellBasedSimulationWithPdes<2> simulator(cell_population, force_collection, pde_and_bc_collection);
        simulator.SetOutputDirectory("TestPostSolveMethod");
        simulator.SetEndTime(2.0/120.0);

        // Set up cell killer and pass into simulation
        OxygenBasedCellKiller<2> killer(&cell_population);
        simulator.AddCellKiller(&killer);

        // Run cell-based simulation
        simulator.Solve();

        // Check the correct solution was obtained
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double radius = norm_2(cell_population.GetLocationOfCellCentre(*cell_iter));
            double analytic_solution = 1 - 0.25*(1 - pow(radius,2.0));

            // Get cell model
            AbstractCellCycleModel* p_abstract_model = cell_iter->GetCellCycleModel();
            SimpleOxygenBasedCellCycleModel* p_oxygen_model = static_cast<SimpleOxygenBasedCellCycleModel*> (p_abstract_model);

            // First part of test - check that PDE solver is working correctly
            TS_ASSERT_DELTA(p_data->GetValue(*cell_iter), analytic_solution, 1e-2);

            // Second part of test - check that each cell's hypoxic duration is correctly updated
            if (p_data->GetValue(*cell_iter) >= p_oxygen_model->GetHypoxicConcentration())
            {
                TS_ASSERT_DELTA(p_oxygen_model->GetCurrentHypoxicDuration(), 0.0, 1e-5);
            }
            else
            {
                TS_ASSERT_DELTA(p_oxygen_model->GetCurrentHypoxicDuration(), 2/120.0, 1e-5);
            }
        }

        // Tidy up
        CellwiseData<2>::Destroy();
    }


    void TestWithOxygen() throw(Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        // Set up mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(STEM);

            //Use non default G1Durations
            p_model->SetStemCellG1Duration(8.0);
            p_model->SetTransitCellG1Duration(8.0);

            CellPtr p_cell(new Cell(p_state, p_model));

            double birth_time = -1.0 - ((double) i/p_mesh->GetNumNodes())*18.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Set up cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetOutputCellPopulationVolumes(true); // record the spheroid radius and apoptotic radius

        // Set up CellwiseData and associate it with the cell population
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);

        // Since values are first passed in to CellwiseData before it is updated in PostSolve(),
        // we need to pass it some initial conditions to avoid memory errors
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, p_mesh->GetNode(i)->GetIndex());
        }

        // Set up PDE
        SimpleUniformSourcePde<2> pde(-0.1);
        double boundary_value = 1.0;
        bool is_neumann_bc = false;
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, boundary_value, is_neumann_bc);
        std::vector<PdeAndBoundaryConditions<2>*> pde_and_bc_collection;
        pde_and_bc_collection.push_back(&pde_and_bc);

        // Set up force law
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.SetCutOffLength(1.5);
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Set up cell-based simulation
        boost::shared_ptr<CellBasedSimulation<2> > p_simulator(new CellBasedSimulationWithPdes<2>(cell_population, force_collection, pde_and_bc_collection));
        TS_ASSERT_EQUALS(p_simulator->GetIdentifier(), "CellBasedSimulationWithPdes-2");
        p_simulator->SetOutputDirectory("CellBasedSimulationWithOxygen");
        p_simulator->SetEndTime(0.5);

        // Set up cell killer and pass into simulation
        OxygenBasedCellKiller<2> killer(&cell_population);
        p_simulator->AddCellKiller(&killer);

        // Run cell-based simulation
        TS_ASSERT_THROWS_NOTHING(p_simulator->Solve());

        // Tidy up
        CellwiseData<2>::Destroy();
    }


    /*
     * This test compares the visualizer output from the previous test
     * with a known file.
     *
     * Note: if the previous test is changed we need to update the file
     * this test refers to.
     */
    void TestVisualizerOutput() throw (Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        // Work out where one of the previous tests wrote its files
        OutputFileHandler handler("CellBasedSimulationWithOxygen", false);
        std::string results_dir = handler.GetOutputDirectoryFullPath() + "results_from_time_0";

        NumericFileComparison comp_nut(results_dir + "/results.vizpdesolution", "cell_based/test/data/CellBasedSimulationWithOxygen/results.vizpdesolution");
        TS_ASSERT(comp_nut.CompareFiles());
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "/results.vizpdesolution cell_based/test/data/CellBasedSimulationWithOxygen/results.vizpdesolution").c_str()), 0);

        NumericFileComparison comp_ele(results_dir + "/results.vizelements", "cell_based/test/data/CellBasedSimulationWithOxygen/results.vizelements");
        TS_ASSERT(comp_ele.CompareFiles());
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "/results.vizelements cell_based/test/data/CellBasedSimulationWithOxygen/results.vizelements").c_str()), 0);

        NumericFileComparison comp_nodes(results_dir + "/results.viznodes", "cell_based/test/data/CellBasedSimulationWithOxygen/results.viznodes");
        TS_ASSERT(comp_nodes.CompareFiles(1e-15));

        NumericFileComparison comp_celltypes(results_dir + "/results.vizcelltypes", "cell_based/test/data/CellBasedSimulationWithOxygen/results.vizcelltypes");
        TS_ASSERT(comp_celltypes.CompareFiles(1e-15));

        TS_ASSERT_EQUALS(system(("diff " + results_dir + "/results.vizsetup cell_based/test/data/CellBasedSimulationWithOxygen/results.vizsetup").c_str()), 0);
    }


    void TestWithPointwiseSource() throw(Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        // Set up mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_state(new ApoptoticCellProperty);

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(STEM);
            CellPtr p_cell(new Cell(p_state, p_model));

            //Use non default G1Durations
            p_model->SetStemCellG1Duration(8.0);
            p_model->SetTransitCellG1Duration(8.0);

            double birth_time = -1.0 - ((double) i/p_mesh->GetNumNodes())*18.0;
            p_cell->SetBirthTime(birth_time);

            // Make the cell apoptotic if near the centre
            double x = p_mesh->GetNode(i)->rGetLocation()[0];
            double y = p_mesh->GetNode(i)->rGetLocation()[1];
            double dist_from_centre = sqrt( (x-2.5)*(x-2.5) + (y-2.5)*(y-2.5) );

            if (dist_from_centre < 1.5)
            {
                p_cell->AddCellProperty(p_apoptotic_state);
            }

            cells.push_back(p_cell);
        }

        // Set up cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetOutputCellPopulationVolumes(true); // record the spheroid radius and apoptotic radius

        // Set up CellwiseData and associate it with the cell population
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);

        // Since values are first passed in to CellwiseData before it is updated in PostSolve(),
        // we need to pass it some initial conditions to avoid memory errors
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, p_mesh->GetNode(i)->GetIndex());
        }

        // Set up PDE
        CellwiseSourcePde<2> pde(cell_population, -0.1);
        double boundary_value = 1.0;
        bool is_neumann_bc = false;
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, boundary_value, is_neumann_bc);
        std::vector<PdeAndBoundaryConditions<2>*> pde_and_bc_collection;
        pde_and_bc_collection.push_back(&pde_and_bc);

        // Set up force law
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.SetCutOffLength(1.5);
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Set up cell-based simulation
        CellBasedSimulationWithPdes<2> simulator(cell_population, force_collection, pde_and_bc_collection);
        simulator.SetOutputDirectory("CellBasedSimulationWithPdes");
        simulator.SetEndTime(0.5);

        // Set up cell killer and pass into simulation
        OxygenBasedCellKiller<2> killer(&cell_population);
        simulator.AddCellKiller(&killer);

        // Run cell-based simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // A few hardcoded tests to check nothing has changed
        std::vector<double> node_5_location = simulator.GetNodeLocation(5);
        TS_ASSERT_DELTA(node_5_location[0], 0.6576, 1e-4);
        TS_ASSERT_DELTA(node_5_location[1], 1.1358, 1e-4);
        TS_ASSERT_DELTA(p_data->GetValue(simulator.rGetCellPopulation().GetCellUsingLocationIndex(5)), 0.9702, 1e-4);

        // Tidy up
        CellwiseData<2>::Destroy();
    }

    void TestWithPointwiseTwoSource() throw(Exception)
    {
		EXIT_IF_PARALLEL; // defined in PetscTools

		// Set up mesh
		unsigned num_cells_depth = 5;
		unsigned num_cells_width = 5;
		HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
		MutableMesh<2,2>* p_mesh = generator.GetMesh();

		// Set up cells
		std::vector<CellPtr> cells;
		boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
		boost::shared_ptr<AbstractCellProperty> p_apoptotic_state(new ApoptoticCellProperty);

		for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
		{
			SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();
			p_model->SetDimension(2);
			p_model->SetCellProliferativeType(STEM);

			//Use non default G1Durations
            p_model->SetStemCellG1Duration(8.0);
	        p_model->SetTransitCellG1Duration(8.0);

			CellPtr p_cell(new Cell(p_state, p_model));
			double birth_time = -1.0 - ((double) i/p_mesh->GetNumNodes())*18.0;
			p_cell->SetBirthTime(birth_time);

			// Make the cell apoptotic if near the centre
			double x = p_mesh->GetNode(i)->rGetLocation()[0];
			double y = p_mesh->GetNode(i)->rGetLocation()[1];
			double dist_from_centre = sqrt( (x-2.5)*(x-2.5) + (y-2.5)*(y-2.5) );

			if (dist_from_centre < 1.5)
			{
				p_cell->AddCellProperty(p_apoptotic_state);
			}

			cells.push_back(p_cell);
		}

		// Set up cell population
		MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
		cell_population.SetOutputCellPopulationVolumes(true); // record the spheroid radius and apoptotic radius

		// Set up CellwiseData and associate it with the cell population
		CellwiseData<2>* p_data = CellwiseData<2>::Instance();
		p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 2);
		p_data->SetCellPopulation(&cell_population);

		// Since values are first passed in to CellwiseData before it is updated in PostSolve(),
		// we need to pass it some initial conditions to avoid memory errors
		for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
		{
			p_data->SetValue(1.0, p_mesh->GetNode(i)->GetIndex(),0);
			p_data->SetValue(1.0, p_mesh->GetNode(i)->GetIndex(),1);
		}

		// Set up PDE
		CellwiseSourcePde<2> pde(cell_population, -0.1);
		double boundary_value = 1.0;
		bool is_neumann_bc = false;
		PdeAndBoundaryConditions<2> pde_and_bc(&pde, boundary_value, is_neumann_bc);
		std::vector<PdeAndBoundaryConditions<2>*> pde_and_bc_collection;
		pde_and_bc_collection.push_back(&pde_and_bc);
		// Set up second PDE
		CellwiseSourcePde<2> pde2(cell_population, -0.8);
		double boundary_value2 = 0.0; //zero flux BC
		bool is_neumann_bc2 = true;
		PdeAndBoundaryConditions<2> pde_and_bc2(&pde2, boundary_value2, is_neumann_bc2);
		pde_and_bc_collection.push_back(&pde_and_bc2);

		// Set up force law
		GeneralisedLinearSpringForce<2> linear_force;
		linear_force.SetCutOffLength(1.5);
		std::vector<AbstractForce<2>*> force_collection;
		force_collection.push_back(&linear_force);

		// Set up cell-based simulation
		CellBasedSimulationWithPdes<2> simulator(cell_population, force_collection, pde_and_bc_collection);
		simulator.SetOutputDirectory("CellBasedSimulationWithPointwiseSource");
		simulator.SetEndTime(0.5);

		// Set up cell killer and pass into simulation
		OxygenBasedCellKiller<2> killer(&cell_population);
		simulator.AddCellKiller(&killer);

		// Run cell-based simulation
		TS_ASSERT_THROWS_NOTHING(simulator.Solve());

		// A few hardcoded tests to check nothing has changed
		std::vector<double> node_5_location = simulator.GetNodeLocation(5);
		TS_ASSERT_DELTA(node_5_location[0], 0.6576, 1e-4);
		TS_ASSERT_DELTA(node_5_location[1], 1.1358, 1e-4);
		TS_ASSERT_DELTA(p_data->GetValue(simulator.rGetCellPopulation().GetCellUsingLocationIndex(5),0), 0.9702, 1e-4);
		TS_ASSERT_DELTA(p_data->GetValue(simulator.rGetCellPopulation().GetCellUsingLocationIndex(5),1), 0.0000, 1e-4);
		TS_ASSERT_LESS_THAN(p_data->GetValue(simulator.rGetCellPopulation().GetCellUsingLocationIndex(5),1),
							p_data->GetValue(simulator.rGetCellPopulation().GetCellUsingLocationIndex(5),0));

		// Tidy up
		CellwiseData<2>::Destroy();
	}

    void TestSpheroidStatistics() throw (Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        // Set up mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_state(new ApoptoticCellProperty);

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(STEM);

            //Use non default G1Durations
            p_model->SetStemCellG1Duration(8.0);
            p_model->SetTransitCellG1Duration(8.0);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetBirthTime(-0.1);

            // Label three neighbouring cells as apoptotic
            if (i==12 || i==13 || i==17)
            {
                p_cell->AddCellProperty(p_apoptotic_state);
            }
            cells.push_back(p_cell);
        }

        // Set up cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetOutputCellPopulationVolumes(true); // record the spheroid radius and apoptotic radius

        // Set up CellwiseData and associate it with the cell population
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, p_mesh->GetNode(i)->GetIndex());
        }

        // Set up PDE
        SimpleUniformSourcePde<2> pde(-0.1);
        double boundary_value = 1.0;
        bool is_neumann_bc = false;
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, boundary_value, is_neumann_bc);
        std::vector<PdeAndBoundaryConditions<2>*> pde_and_bc_collection;
        pde_and_bc_collection.push_back(&pde_and_bc);

        // Set up force law
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.SetCutOffLength(1.5);
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Set up cell-based simulation
        CellBasedSimulationWithPdes<2> simulator(cell_population, force_collection, pde_and_bc_collection);
        simulator.SetOutputDirectory("TestSpheroidStatistics");
        simulator.SetEndTime(1.0/120.0);
        simulator.SetWriteAverageRadialPdeSolution(5);

        // Add an oxygen-dependent cell killer to the cell-based simulation
        OxygenBasedCellKiller<2> killer(&cell_population);
        simulator.AddCellKiller(&killer);

        // Run the cell-based simulation for one timestep
        simulator.Solve();

        // Just check that we do indeed have three apoptotic cells
        unsigned num_apoptotic_cells = 0;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            if (cell_iter->HasCellProperty<ApoptoticCellProperty>())
            {
                num_apoptotic_cells++;
            }
        }
        TS_ASSERT_EQUALS(num_apoptotic_cells, 3u);

        /**
         * We have 25 cells. Adding up the boundary cell areas, we
         * should have the equivalent area of 16 full regular hexagonal
         * cells.
         *
         * The area of a single hexagonal cell is sqrt(3)/2, so the
         * correct spheroid radius is given by sqrt((16*sqrt(3)/2)/pi).
         *
         * Since there are 3 apoptotic cells, the correct apoptotic radius
         * is given by sqrt((3*sqrt(3)/2)/pi).
         */

        // Work out where the previous test wrote its files
        OutputFileHandler handler("TestSpheroidStatistics", false);
        std::string areas_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/cellpopulationareas.dat";
        TS_ASSERT_EQUALS(system(("diff " + areas_results_file + " cell_based/test/data/TestSpheroidStatistics/cellpopulationareas.dat").c_str()), 0);

        std::string dist_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/radial_dist.dat";
        TS_ASSERT_EQUALS(system(("diff " + dist_results_file + " cell_based/test/data/TestSpheroidStatistics/radial_dist.dat").c_str()), 0);

        // Coverage
        TS_ASSERT_THROWS_NOTHING(simulator.WriteAverageRadialPdeSolution(SimulationTime::Instance()->GetTime(),5));

        // Tidy up
        CellwiseData<2>::Destroy();
    }


    void TestCoarseSourceMesh() throw(Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        // Set up mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(STEM);

            //Use non default G1Durations
			p_model->SetStemCellG1Duration(8.0);
			p_model->SetTransitCellG1Duration(8.0);


            CellPtr p_cell(new Cell(p_state, p_model));
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*18.0;
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Set up cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up CellwiseData and associate it with the cell population
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 2);
        p_data->SetCellPopulation(&cell_population);

        // Since values are first passed in to CellwiseData before it is updated in PostSolve(),
        // we need to pass it some initial conditions to avoid memory errors
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, p_mesh->GetNode(i)->GetIndex(),0);
            p_data->SetValue(1.0, p_mesh->GetNode(i)->GetIndex(),1);
        }

        // Set up PDE
        AveragedSourcePde<2> pde(cell_population, -0.1);
        double boundary_value = 1.0;
        bool is_neumann_bc = false;
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, boundary_value, is_neumann_bc);

        // Set up second PDE
        AveragedSourcePde<2> pde2(cell_population, -0.5);
        double boundary_value2 = 1.0;
        bool is_neumann_bc2 = false;
        PdeAndBoundaryConditions<2> pde_and_bc2(&pde2, boundary_value2, is_neumann_bc2);

        std::vector<PdeAndBoundaryConditions<2>*> pde_and_bc_collection;
        pde_and_bc_collection.push_back(&pde_and_bc);
        pde_and_bc_collection.push_back(&pde_and_bc2);

        // Set up force law
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.SetCutOffLength(1.5);
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Set up cell-based simulation
        CellBasedSimulationWithPdes<2> simulator(cell_population, force_collection, pde_and_bc_collection);
        simulator.SetOutputDirectory("TestCoarseSourceMesh");
        simulator.SetEndTime(0.05);

        // Coverage
        simulator.SetPdeAndBcCollection(pde_and_bc_collection);

        // Set up cell killer and pass into simulation
        OxygenBasedCellKiller<2> killer(&cell_population);
        simulator.AddCellKiller(&killer);

        // Test creation of mpCoarsePdeMesh
        simulator.UseCoarsePdeMesh(10.0);

        // Find centre of cell population
        c_vector<double,2> centre_of_cell_population = zero_vector<double>(2);

        for (unsigned i=0; i<simulator.rGetCellPopulation().GetNumNodes(); i++)
        {
            centre_of_cell_population += simulator.rGetCellPopulation().GetNode(i)->rGetLocation();
        }
        centre_of_cell_population /= simulator.rGetCellPopulation().GetNumNodes();

        // Find centre of coarse PDE mesh
        c_vector<double,2> centre_of_coarse_pde_mesh = zero_vector<double>(2);

        for (unsigned i=0; i<simulator.mpCoarsePdeMesh->GetNumNodes(); i++)
        {
            centre_of_coarse_pde_mesh += simulator.mpCoarsePdeMesh->GetNode(i)->rGetLocation();
        }
        centre_of_coarse_pde_mesh /= simulator.mpCoarsePdeMesh->GetNumNodes();

        // Test that the two centres match
        TS_ASSERT_DELTA(centre_of_cell_population[0], centre_of_coarse_pde_mesh[0], 1e-4);
        TS_ASSERT_DELTA(centre_of_cell_population[1], centre_of_coarse_pde_mesh[1], 1e-4);

        // Test FindCoarseElementContainingCell and initialisation of mCellPdeElementMap

        simulator.InitialiseCoarsePdeMesh(); // coverage

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
            cell_iter != cell_population.End();
            ++cell_iter)
        {
            unsigned containing_element_index = simulator.mCellPdeElementMap[*cell_iter];
            TS_ASSERT_LESS_THAN(containing_element_index, simulator.mpCoarsePdeMesh->GetNumElements());
            TS_ASSERT_EQUALS(containing_element_index, simulator.FindCoarseElementContainingCell(*cell_iter));
        }

        // Run cell-based simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
        TS_ASSERT(simulator.mpCoarsePdeMesh != NULL);

        ReplicatableVector pde_solution0(simulator.GetCurrentPdeSolution(0));
        ReplicatableVector pde_solution1(simulator.GetCurrentPdeSolution(1));

        TS_ASSERT_EQUALS(pde_solution0.GetSize(),pde_solution1.GetSize());

        // Test the nutrient concentration is equal to 1.0 at each coarse mesh node far from the cells
        for (unsigned i=0; i<pde_solution0.GetSize(); i++)
        {
            c_vector<double,2> centre;
            centre(0) = 2.5; // assuming 5 by 5 honeycomb mesh
            centre(1) = 2.5;
            c_vector<double,2> posn = simulator.mpCoarsePdeMesh->GetNode(i)->rGetLocation();
            double dist = norm_2(centre - posn);
            double u0 = pde_solution0[i];
            double u1 = pde_solution1[i];

            if (dist > 4.0)
            {
                TS_ASSERT_DELTA(u0, 1.0, 1e-5);
                TS_ASSERT_DELTA(u1, 1.0, 1e-5);
            }
        }

        /*
         * Loop over cells, find the coarse mesh element containing it, then
         * check the interpolated PDE solution is between the min and max of
         * the PDE solution on the nodes of that element.
         */
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
            cell_iter != cell_population.End();
            ++cell_iter)
        {
            unsigned elem_index = simulator.mpCoarsePdeMesh->GetContainingElementIndex(cell_population.GetLocationOfCellCentre(*cell_iter));
            Element<2,2>* p_element = simulator.mpCoarsePdeMesh->GetElement(elem_index);


            double max0 = std::max(pde_solution0[p_element->GetNodeGlobalIndex(0)],
                                  pde_solution0[p_element->GetNodeGlobalIndex(1)]);
            max0 = std::max(max0, pde_solution0[p_element->GetNodeGlobalIndex(2)]);

            double max1 = std::max(pde_solution1[p_element->GetNodeGlobalIndex(0)],
                                   pde_solution1[p_element->GetNodeGlobalIndex(1)]);
            max1 = std::max(max1, pde_solution1[p_element->GetNodeGlobalIndex(2)]);

            double min0 = std::min(pde_solution0[p_element->GetNodeGlobalIndex(0)],
                                  pde_solution0[p_element->GetNodeGlobalIndex(1)]);
            min0 = std::min(min0, pde_solution0[p_element->GetNodeGlobalIndex(2)]);

            double min1 = std::min(pde_solution1[p_element->GetNodeGlobalIndex(0)],
                                  pde_solution1[p_element->GetNodeGlobalIndex(1)]);
            min1 = std::min(min1, pde_solution1[p_element->GetNodeGlobalIndex(2)]);

            double value0_at_cell = CellwiseData<2>::Instance()->GetValue(*cell_iter, 0);
            double value1_at_cell = CellwiseData<2>::Instance()->GetValue(*cell_iter, 1);

            TS_ASSERT_LESS_THAN_EQUALS(value1_at_cell, value0_at_cell);
            TS_ASSERT_LESS_THAN_EQUALS(min0, value0_at_cell);
            TS_ASSERT_LESS_THAN_EQUALS(value0_at_cell, max0);
            TS_ASSERT_LESS_THAN_EQUALS(min1, value1_at_cell);
            TS_ASSERT_LESS_THAN_EQUALS(value1_at_cell, max1);
        }

        // Tidy up
        CellwiseData<2>::Destroy();
    }


    void TestCoarseSourceMeshWithNeumannIsNotImplemented() throw(Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        // Set up mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(STEM);

            //Use non default G1Durations
			p_model->SetStemCellG1Duration(8.0);
			p_model->SetTransitCellG1Duration(8.0);

            CellPtr p_cell(new Cell(p_state, p_model));
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*18.0;
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Set up cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up CellwiseData and associate it with the cell population
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 2);
        p_data->SetCellPopulation(&cell_population);

        // Since values are first passed in to CellwiseData before it is updated in PostSolve(),
        // we need to pass it some initial conditions to avoid memory errors
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, p_mesh->GetNode(i)->GetIndex(),0);
            p_data->SetValue(1.0, p_mesh->GetNode(i)->GetIndex(),1);
        }

        // Set up PDE
        AveragedSourcePde<2> pde(cell_population, -0.1);
        double boundary_value = 0.0;
        bool is_neumann_bc = true;
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, boundary_value, is_neumann_bc);

        // Set up second PDE
        AveragedSourcePde<2> pde2(cell_population, -0.5);
        double boundary_value2 = 0.0;
        bool is_neumann_bc2 = true;
        PdeAndBoundaryConditions<2> pde_and_bc2(&pde2, boundary_value2, is_neumann_bc2);

        std::vector<PdeAndBoundaryConditions<2>*> pde_and_bc_collection;
        pde_and_bc_collection.push_back(&pde_and_bc);
        pde_and_bc_collection.push_back(&pde_and_bc2);

        // Set up force law
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.SetCutOffLength(1.5);
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Set up cell-based simulation
        CellBasedSimulationWithPdes<2> simulator(cell_population, force_collection, pde_and_bc_collection);
        simulator.SetOutputDirectory("TestCoarseSourceMesh");
        simulator.SetEndTime(0.05);

        // Coverage
        simulator.SetPdeAndBcCollection(pde_and_bc_collection);

        // Set up cell killer and pass into simulation
        OxygenBasedCellKiller<2> killer(&cell_population);
        simulator.AddCellKiller(&killer);

        // Test creation of mpCoarsePdeMesh
        simulator.UseCoarsePdeMesh(10.0);

        // Find centre of cell population
        c_vector<double,2> centre_of_cell_population = zero_vector<double>(2);

        for (unsigned i=0; i<simulator.rGetCellPopulation().GetNumNodes(); i++)
        {
            centre_of_cell_population += simulator.rGetCellPopulation().GetNode(i)->rGetLocation();
        }
        centre_of_cell_population /= simulator.rGetCellPopulation().GetNumNodes();

        // Find centre of coarse PDE mesh
        c_vector<double,2> centre_of_coarse_pde_mesh = zero_vector<double>(2);

        for (unsigned i=0; i<simulator.mpCoarsePdeMesh->GetNumNodes(); i++)
        {
            centre_of_coarse_pde_mesh += simulator.mpCoarsePdeMesh->GetNode(i)->rGetLocation();
        }
        centre_of_coarse_pde_mesh /= simulator.mpCoarsePdeMesh->GetNumNodes();

        // Test that the two centres match
        TS_ASSERT_DELTA(centre_of_cell_population[0], centre_of_coarse_pde_mesh[0], 1e-4);
        TS_ASSERT_DELTA(centre_of_cell_population[1], centre_of_coarse_pde_mesh[1], 1e-4);

        // Test FindCoarseElementContainingCell and initialisation of mCellPdeElementMap

        simulator.InitialiseCoarsePdeMesh(); // coverage

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
            cell_iter != cell_population.End();
            ++cell_iter)
        {
            unsigned containing_element_index = simulator.mCellPdeElementMap[*cell_iter];
            TS_ASSERT_LESS_THAN(containing_element_index, simulator.mpCoarsePdeMesh->GetNumElements());
            TS_ASSERT_EQUALS(containing_element_index, simulator.FindCoarseElementContainingCell(*cell_iter));
        }

        // Run cell-based simulation
        TS_ASSERT_THROWS_THIS(simulator.Solve(), "Neumann BCs not yet implemented when using a coarse PDE mesh");



        // Tidy up
        CellwiseData<2>::Destroy();
    }

    void TestCoarseSourceMeshBoundaryConditionImplementation() throw(Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        // Create a cigar-shaped mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        MutableMesh<2,2>* p_mesh = new MutableMesh<2,2>;
        p_mesh->ConstructFromMeshReader(mesh_reader);
        p_mesh->Scale(5.0,1.0);

        // Set up cells
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(STEM);

            //Use non default G1Durations
			p_model->SetStemCellG1Duration(8.0);
			p_model->SetTransitCellG1Duration(8.0);


            CellPtr p_cell(new Cell(p_state, p_model));
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*18.0;
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Set up cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up CellwiseData and associate it with the cell population
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);

        // Since values are first passed in to CellwiseData before it is updated in PostSolve(),
        // we need to pass it some initial conditions to avoid memory errors
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, p_mesh->GetNode(i)->GetIndex());
        }

        // Set up PDE
        AveragedSourcePde<2> pde(cell_population, -0.01);
        double boundary_value = 1.0;
        bool is_neumann_bc = false;
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, boundary_value, is_neumann_bc);
        std::vector<PdeAndBoundaryConditions<2>*> pde_and_bc_collection;
        pde_and_bc_collection.push_back(&pde_and_bc);

        // Set up force law
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.SetCutOffLength(1.5);
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Set up cell-based simulation to use a coarse PDE mesh
        CellBasedSimulationWithPdes<2> simulator(cell_population, force_collection, pde_and_bc_collection);
        simulator.SetOutputDirectory("TestCoarseNutrientMeshBoundaryConditionImplementation");
        simulator.SetEndTime(0.01);
        simulator.UseCoarsePdeMesh(2.0);

        // Run cell-based simulation
        simulator.Solve();

        // Test that boundary cells experience the right boundary condition
        for (AbstractCellPopulation<2>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
             cell_iter != simulator.rGetCellPopulation().End();
             ++cell_iter)
        {
            if ( (static_cast<AbstractCentreBasedCellPopulation<2>*>(&(simulator.rGetCellPopulation())))->GetNodeCorrespondingToCell(*cell_iter)->IsBoundaryNode() )
            {
                TS_ASSERT_DELTA(p_data->GetValue(*cell_iter), 1.0, 1e-1);
            }
        }

        // Tidy up
        CellwiseData<2>::Destroy();
        delete p_mesh;
    }


    void TestArchivingWithSimplePde() throw (Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        // Set up mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(STEM);

            //Use non default G1Durations
			p_model->SetStemCellG1Duration(8.0);
			p_model->SetTransitCellG1Duration(8.0);


            CellPtr p_cell(new Cell(p_state, p_model));
            double birth_time = -1.0 - ((double) i/p_mesh->GetNumNodes())*18.0;
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Set up cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up CellwiseData and associate it with the cell population
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);

        // Since values are first passed in to CellwiseData before it is updated in PostSolve(),
        // we need to pass it some initial conditions to avoid memory errors
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, p_mesh->GetNode(i)->GetIndex());
        }

        // Set up PDE
        SimpleUniformSourcePde<2> pde(-0.1);
        double boundary_value = 1.0;
        bool is_neumann_bc = false;
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, boundary_value, is_neumann_bc);
        std::vector<PdeAndBoundaryConditions<2>*> pde_and_bc_collection;
        pde_and_bc_collection.push_back(&pde_and_bc);

        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.SetCutOffLength(1.5);
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Set up cell-based simulation
        CellBasedSimulationWithPdes<2> simulator(cell_population, force_collection, pde_and_bc_collection);
        simulator.SetOutputDirectory("CellBasedSimulationWithPdesSaveAndLoad");
        simulator.SetEndTime(0.2);

        // Set up cell killer and pass into simulation
        OxygenBasedCellKiller<2> killer(&cell_population);
        simulator.AddCellKiller(&killer);

        // Run cell-based simulation
        simulator.Solve();

        // Save cell-based simulation
        CellBasedSimulationArchiver<2, CellBasedSimulationWithPdes<2> >::Save(&simulator);

        CellBasedSimulationWithPdes<2>* p_simulator
            = CellBasedSimulationArchiver<2, CellBasedSimulationWithPdes<2> >::Load("CellBasedSimulationWithPdesSaveAndLoad", 0.2);

        p_simulator->SetPdeAndBcCollection(pde_and_bc_collection);
        p_simulator->SetEndTime(0.5);
        p_simulator->Solve();

        // These results are from time 0.5 in TestWithOxygen.
        std::vector<double> node_5_location = p_simulator->GetNodeLocation(5);
        TS_ASSERT_DELTA(node_5_location[0], 0.4968, 1e-4);
        TS_ASSERT_DELTA(node_5_location[1], 0.8635, 1e-4);

        std::vector<double> node_15_location = p_simulator->GetNodeLocation(15);
        TS_ASSERT_DELTA(node_15_location[0], 0.4976, 1e-4);
        TS_ASSERT_DELTA(node_15_location[1], 2.5977, 1e-4);

        // Test CellwiseData was set up correctly
        TS_ASSERT_EQUALS(CellwiseData<2>::Instance()->IsSetUp(), true);

        // Test the CellwiseData result
        TS_ASSERT_DELTA(p_data->GetValue(p_simulator->rGetCellPopulation().GetCellUsingLocationIndex(5)), 0.9604, 1e-4);
        TS_ASSERT_DELTA(p_data->GetValue(p_simulator->rGetCellPopulation().GetCellUsingLocationIndex(15)), 0.9584, 1e-4);

        // Run cell-based simulation
        delete p_simulator;
        CellwiseData<2>::Destroy();
    }


    /**
     * This test demonstrates how to archive a CellBasedSimulationWithPdes
     * in the case where the PDE has the cell population as a member variable.
     */
    void TestArchivingWithCellwisePde() throw (Exception)
    {
        if (!PetscTools::IsSequential())
        {
            TS_TRACE("This test does not pass in parallel yet.");
            return;
        }

        std::string output_directory = "TestArchivingWithCellwisePde";
        double end_time = 0.1;

        // Set up mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(STEM);

            //Use non default G1Durations
			p_model->SetStemCellG1Duration(8.0);
			p_model->SetTransitCellG1Duration(8.0);

            CellPtr p_cell(new Cell(p_state, p_model));
            double birth_time = -1.0 - ((double) i/p_mesh->GetNumNodes())*18.0;
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Set up cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up CellwiseData and associate it with the cell population
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);

        // Set initial conditions for CellwiseData (needed to avoid memory errors)
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, p_mesh->GetNode(i)->GetIndex());
        }

        // Set up PDE
        CellwiseSourcePde<2> pde(cell_population, -0.03);
        double boundary_value = 1.0;
        bool is_neumann_bc = false;
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, boundary_value, is_neumann_bc);
        std::vector<PdeAndBoundaryConditions<2>*> pde_and_bc_collection;
        pde_and_bc_collection.push_back(&pde_and_bc);

        // Set up mechanics system
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.SetCutOffLength(3.0);
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        // Set up cell-based simulation
        CellBasedSimulationWithPdes<2> simulator(cell_population, force_collection, pde_and_bc_collection);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetEndTime(end_time);

        // Run cell-based simulation
        simulator.Solve();

        // Save cell-based simulation
        CellBasedSimulationArchiver<2, CellBasedSimulationWithPdes<2> >::Save(&simulator);

        // Load simulation
        CellBasedSimulationWithPdes<2>* p_simulator
            = CellBasedSimulationArchiver<2, CellBasedSimulationWithPdes<2> >::Load(output_directory, end_time);

        /**
         * In this case, the PDE had a reference to the cell population. To avoid a
         * segmentation fault, we need to first get the archived cell_population, pass
         * it in to a new instance of our PDE, taking care to use the same
         * consumption rate as before. We then pass this PDE into the cell-based
         * simulation.
         */
        MeshBasedCellPopulation<2>* p_cell_population = static_cast<MeshBasedCellPopulation<2>*>(&(p_simulator->rGetCellPopulation()));
        CellwiseSourcePde<2> pde2(*p_cell_population, -0.03);
        double boundary_value2 = 1.0;
        bool is_neumann_bc2 = false;
        PdeAndBoundaryConditions<2> pde_and_bc2(&pde2, boundary_value2, is_neumann_bc2);
        std::vector<PdeAndBoundaryConditions<2>*> pde_and_bc_collection2;
        pde_and_bc_collection2.push_back(&pde_and_bc2);

        p_simulator->SetPdeAndBcCollection(pde_and_bc_collection2);
        p_simulator->SetEndTime(2.0*end_time);

        // Run cell-based simulation
        TS_ASSERT_THROWS_NOTHING(p_simulator->Solve());

        // Tidy up
        delete p_simulator;
        CellwiseData<2>::Destroy();
    }


    void Test3DCellBasedSimulationWithPdes() throw(Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        MutableMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TrianglesMeshWriter<3,3> mesh_writer("TestSolveMethodSpheroidSimulation3DMesh", "StartMesh");
        mesh_writer.WriteFilesUsingMesh(mesh);

        // Set up cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> generator;
        generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Set some model parameters for the cell cycle model
		for(unsigned index=0; index < cells.size(); index++)
		{
			cells[index]->GetCellCycleModel()->SetTransitCellG1Duration(8.0);
			cells[index]->GetCellCycleModel()->SetStemCellG1Duration(8.0);
		}

        // Set up cell population
        MeshBasedCellPopulation<3> cell_population(mesh, cells);

        // Set up CellwiseData and associate it with the cell population
        CellwiseData<3>* p_data = CellwiseData<3>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);

        // Since values are first passed in to CellwiseData before it is updated in PostSolve(),
        // we need to pass it some initial conditions to avoid memory errors
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, mesh.GetNode(i)->GetIndex());
        }

        // Set up PDE
        SimpleUniformSourcePde<3> pde(-0.03);
        double boundary_value = 1.0;
        bool is_neumann_bc = false;
        PdeAndBoundaryConditions<3> pde_and_bc(&pde, boundary_value, is_neumann_bc);
        std::vector<PdeAndBoundaryConditions<3>*> pde_and_bc_collection;
        pde_and_bc_collection.push_back(&pde_and_bc);

        // Set up force law
        GeneralisedLinearSpringForce<3> linear_force;
        linear_force.SetCutOffLength(1.5);
        std::vector<AbstractForce<3>*> force_collection;
        force_collection.push_back(&linear_force);

        // Set up cell-based simulation
        CellBasedSimulationWithPdes<3> simulator(cell_population, force_collection, pde_and_bc_collection);
        simulator.SetOutputDirectory("CellBasedSimulationWithOxygen3d");
        simulator.SetEndTime(0.5);

        // Set up cell killer and pass into simulation
        AbstractCellKiller<3>* p_killer = new OxygenBasedCellKiller<3>(&cell_population);
        simulator.AddCellKiller(p_killer);

        // Coverage
        TS_ASSERT_THROWS_THIS(simulator.CreateCoarsePdeMesh(10.0), "This method is only implemented in 2D");

        // Run cell-based simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // Tidy up
        CellwiseData<3>::Destroy();

        delete p_killer;
    }

    void TestCellBasedSimulationWithPdesParameterOutputMethods() throw (Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

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

        // Set up PDE
        CellwiseSourcePde<2> pde(cell_population, -0.03);
        double boundary_value = 1.0;
        bool is_neumann_bc = false;
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, boundary_value, is_neumann_bc);
        std::vector<PdeAndBoundaryConditions<2>*> pde_and_bc_collection;
        pde_and_bc_collection.push_back(&pde_and_bc);

        // Create a force law Collection
		GeneralisedLinearSpringForce<2> linear_force;
		std::vector<AbstractForce<2>*> force_collection;
		force_collection.push_back(&linear_force);

        // Set up cell-based simulation
		CellBasedSimulationWithPdes<2> simulator(cell_population, force_collection, pde_and_bc_collection);

        TS_ASSERT_EQUALS(simulator.GetIdentifier(), "CellBasedSimulationWithPdes-2");

		std::string output_directory = "TestCellBasedSimulationOutputParameters";
		OutputFileHandler output_file_handler(output_directory, false);
		out_stream parameter_file = output_file_handler.OpenOutputFile("cell_based_sim_with_pde_results.parameters");
		simulator.OutputSimulationParameters(parameter_file);
		parameter_file->close();

		std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();
		TS_ASSERT_EQUALS(system(("diff " + results_dir + "cell_based_sim_with_pde_results.parameters  cell_based/test/data/TestCellBasedSimulationOutputParameters/cell_based_sim_with_pde_results.parameters").c_str()), 0);

		///\todo check output of simulator.OutputSimulationSetup();
    }
};

#endif /*TESTCELLBASEDSIMULATIONWITHPDES_HPP_*/
