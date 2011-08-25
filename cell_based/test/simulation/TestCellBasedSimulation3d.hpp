/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef TESTCELLBASEDSIMULATION3D_HPP_
#define TESTCELLBASEDSIMULATION3D_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "CellsGenerator.hpp"
#include "TrianglesMeshReader.hpp"
#include "CellBasedSimulation.hpp"
#include "TrianglesMeshWriter.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "Debug.hpp"

class TestCellBasedSimulation3d : public AbstractCellBasedTestSuite
{
private:
    double mLocationGhosts;
    double mLocationWithoutGhosts;

    MutableMesh<3,3>* Make3dMesh(unsigned width=3, unsigned height=3, unsigned depth=3)
    {
        MutableMesh<3,3>* p_mesh = new MutableMesh<3,3>;
        p_mesh->ConstructCuboid(width, height, depth);
        TrianglesMeshWriter<3,3> mesh_writer("", "3dSpringMesh");
        mesh_writer.WriteFilesUsingMesh(*p_mesh);

        return p_mesh;
    }

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

    void TestDoCellBirth() throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_1626_elements");
        MutableMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        for (unsigned i=0; i<cells.size(); i++)
        {
            dynamic_cast<FixedDurationGenerationBasedCellCycleModel*>(cells[i]->GetCellCycleModel())->SetGeneration(0);
            cells[i]->SetBirthTime(0.0);
        }
        cells[50]->SetBirthTime(-50.0);

        MeshBasedCellPopulation<3> cell_population(mesh, cells);

        CellBasedSimulation<3> simulator(cell_population);

        unsigned num_births = simulator.DoCellBirth();

        TS_ASSERT_EQUALS(num_births, 1u);
    }

    void TestBirthOccursDuringSolve() throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_Single_tetrahedron_element");

        MutableMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        for (unsigned i=0; i<cells.size(); i++)
        {
            dynamic_cast<FixedDurationGenerationBasedCellCycleModel*>(cells[i]->GetCellCycleModel())->SetGeneration(0);
            cells[i]->SetBirthTime(0.0);
        }
        cells[mesh.GetNumNodes()-1]->SetBirthTime(-50.0);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 1u);

        MeshBasedCellPopulation<3> cell_population(mesh, cells);

        CellBasedSimulation<3> simulator(cell_population);

        TrianglesMeshWriter<3,3> mesh_writer1("Test3DCellBirth", "StartMesh");
        mesh_writer1.WriteFilesUsingMesh(mesh);

        simulator.SetOutputDirectory("Test3DCellBirth");
        simulator.SetEndTime(1.0);

        // Create a force law and pass it to the simulation
        GeneralisedLinearSpringForce<3> linear_force;
        linear_force.SetCutOffLength(1.5);
        simulator.AddForce(&linear_force);

        simulator.Solve();

        // Start with a single 3D tetrahedron, add one node get two tetrahedral elements.
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);

        TrianglesMeshWriter<3,3> mesh_writer2("Test3DCellBirth", "EndMesh", false);
        mesh_writer2.WriteFilesUsingMesh(mesh);
    }

    void TestSolveMethodSpheroidSimulation3D() throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        MutableMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TrianglesMeshWriter<3,3> mesh_writer("TestSolveMethodSpheroidSimulation3DMesh", "StartMesh");
        mesh_writer.WriteFilesUsingMesh(mesh);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        for (unsigned i=0; i<cells.size(); i++)
        {
            dynamic_cast<FixedDurationGenerationBasedCellCycleModel*>(cells[i]->GetCellCycleModel())->SetGeneration(0);
        }

        MeshBasedCellPopulation<3> cell_population(mesh, cells);

        CellBasedSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("TestSolveMethodSpheroidSimulation3D");

        // Create a force law and pass it to the simulation
        GeneralisedLinearSpringForce<3> linear_force;
        linear_force.SetCutOffLength(1.5);
        simulator.AddForce(&linear_force);

        // Test SetSamplingTimestepMultiple method
        TS_ASSERT_EQUALS(simulator.mSamplingTimestepMultiple, 1u);
        simulator.SetSamplingTimestepMultiple(2);
        TS_ASSERT_EQUALS(simulator.mSamplingTimestepMultiple, 2u);

        // Uncommenting this line calls an error in accessing nodes in the vertex elements #
        //cell_population.SetOutputVoronoiData(true);

        simulator.SetEndTime(0.1);
        simulator.Solve();

        TrianglesMeshWriter<3,3> mesh_writer2("TestSolveMethodSpheroidSimulation3DMesh", "EndMesh", false);
        mesh_writer2.WriteFilesUsingMesh(mesh);
    }

    void TestGhostNodesSpheroidSimulation3DandSave() throw (Exception)
    {
        unsigned width = 3;
        unsigned height = 3;
        unsigned depth = 3;

        MutableMesh<3,3>* p_mesh = Make3dMesh(width, height, depth);
        TrianglesMeshWriter<3,3> mesh_writer("TestGhostNodesSpheroidSimulation3D", "StartMesh");
        mesh_writer.WriteFilesUsingMesh(*p_mesh);

        c_vector<double, 3> spheroid_centre;
        spheroid_centre[0] = 0.5*((double) width);
        spheroid_centre[1] = 0.5*((double) height);
        spheroid_centre[2] = 0.5*((double) depth);

        // Set up cells by iterating through the mesh nodes
        unsigned num_nodes = p_mesh->GetNumAllNodes();

        PRINT_VARIABLE(num_nodes);

        std::vector<CellPtr> cells;
        std::vector<CellPtr> cells2;
        std::vector<unsigned> location_indices;

        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<num_nodes; i++)
        {
            c_vector<double, 3> node_location = p_mesh->GetNode(i)->rGetLocation();

            unsigned min_spatial_dimension;
            if (width <= height && width <= depth)
            {
                min_spatial_dimension = width;
            }
            else
            {
                if (height <= depth)
                {
                    min_spatial_dimension = height;
                }
                else
                {
                    min_spatial_dimension = depth;
                }
            }

            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(STEM);
            p_model->SetGeneration(0);
            CellPtr p_cell(new Cell(p_state, p_model));

            p_cell->SetBirthTime(-RandomNumberGenerator::Instance()->ranf()*
                                 (p_model->GetStemCellG1Duration() + p_model->GetSG2MDuration()));

            cells2.push_back(p_cell);

            if (norm_2(node_location - spheroid_centre) <= 0.5*sqrt(3)*1.01*((double) min_spatial_dimension)/3.0)
            {
                location_indices.push_back(i);
                cells.push_back(p_cell);
            }
        }

        TS_ASSERT_EQUALS(location_indices.size(), cells.size());
        TS_ASSERT_LESS_THAN(location_indices.size(), num_nodes);
        TS_ASSERT_EQUALS(location_indices.size(), 8u);

        // Test Save with a MeshBasedCellPopulationWithGhostNodes
        MeshBasedCellPopulationWithGhostNodes<3> cell_population(*p_mesh, cells, location_indices);

        cell_population.SetOutputVoronoiData(true);                           // Output Voronoi data
        cell_population.SetOutputCellProliferativeTypes(true);
        cell_population.SetOutputCellMutationStates(true);
        cell_population.SetOutputCellVolumes(true);

        CellBasedSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("TestGhostNodesSpheroidSimulation3D");
        simulator.SetEndTime(0.1);

        // Create a force law and pass it to the simulation
        GeneralisedLinearSpringForce<3> linear_force;
        linear_force.SetCutOffLength(1.0);
        simulator.AddForce(&linear_force);

        simulator.Solve();
        CellBasedSimulationArchiver<3, CellBasedSimulation<3> >::Save(&simulator);

        // To check consistency with for test below
        mLocationGhosts = p_mesh->GetNode(23)->rGetLocation()[2];

        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);

        // Test Save with a MeshBasedCellPopulation - one cell born during this

        MeshBasedCellPopulationWithGhostNodes<3> cell_population2(*p_mesh, cells2);

        cell_population2.SetOutputVoronoiData(true);                           // Output Voronoi data
        cell_population2.SetOutputCellProliferativeTypes(true);
        cell_population2.SetOutputCellMutationStates(true);
        cell_population2.SetOutputCellVolumes(true);

        CellBasedSimulation<3> simulator2(cell_population2);
        simulator2.SetOutputDirectory("TestGhostNodesSpheroidSimulation3DNoGhosts");
        simulator2.SetEndTime(0.1);

        // Pass force_law to the simulation
        simulator2.AddForce(&linear_force);

        simulator2.Solve();
        CellBasedSimulationArchiver<3, CellBasedSimulation<3> >::Save(&simulator2);

        // To check consistency with for test below
        mLocationWithoutGhosts = p_mesh->GetNode(23)->rGetLocation()[2];
        delete p_mesh;
    }

    void TestLoadOf3DSimulation() throw (Exception)
    {
        {
            // With ghost nodes - 56 ghosts 8 real cells
            CellBasedSimulation<3>* p_simulator = CellBasedSimulationArchiver<3, CellBasedSimulation<3> >::Load("TestGhostNodesSpheroidSimulation3D", 0.1);
            unsigned num_cells = p_simulator->rGetCellPopulation().GetNumRealCells();

            TS_ASSERT_EQUALS(num_cells, 8u);
            TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 0.1, 1e-9);
            TS_ASSERT_DELTA(p_simulator->rGetCellPopulation().GetNode(23)->rGetLocation()[2], mLocationGhosts, 1e-6);

            delete p_simulator;
        }

        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);

        {
            // Without ghost nodes - all 65 are real cells
            CellBasedSimulation<3>* p_simulator = CellBasedSimulationArchiver<3, CellBasedSimulation<3> >::Load("TestGhostNodesSpheroidSimulation3DNoGhosts", 0.1);
            unsigned num_cells = p_simulator->rGetCellPopulation().GetNumRealCells();

            TS_ASSERT_EQUALS(num_cells, 65u);
            TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 0.1, 1e-9);
            CellPtr p_cell = p_simulator->rGetCellPopulation().GetCellUsingLocationIndex(23u);
            TS_ASSERT_DELTA(p_simulator->rGetCellPopulation().GetNode(23)->rGetLocation()[2], mLocationWithoutGhosts, 1e-6);

            delete p_simulator;
        }
    }

    /* A test that fails with the following error:
     *
     * scons: building terminated because of errors.
     * cell_based/build/optimised/src/population/MeshBasedCellPopulationWithGhostNodes.o: In function
     * `MeshBasedCellPopulationWithGhostNodes<1u>::WriteVtkResultsToFile()':
     * MeshBasedCellPopulationWithGhostNodes.cpp:(.text._ZN37MeshBasedCellPopulationWithGhostNodesILj1EE21WriteVtkResultsToFileEv
     * [MeshBasedCellPopulationWithGhostNodes<1u>::WriteVtkResultsToFile()]+0xb8e): undefined reference to `MeshBasedCellPopulation<1u>::
     * TessellateIfNeeded()'
     *
     */
    void xTestSimpleCubeFailingTest() throw (Exception)
    {
        unsigned width = 5;
        unsigned height = 5;
        unsigned depth = 4;

        MutableMesh<3,3>* p_mesh = Make3dMesh(width, height, depth);

        // Set up cells by iterating through the mesh nodes
        unsigned num_nodes = p_mesh->GetNumAllNodes();

        unsigned num_epithelial_cells = (width+1)*(height+1);

    	boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);

    	std::vector<CellPtr> cells;
		for (unsigned i=0; i<num_nodes; i++)
		{
			FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(DIFFERENTIATED);
            p_model->SetDimension(3);

            CellPtr p_differentiated_cell(new Cell(p_state, p_model));
			cells.push_back(p_differentiated_cell);
        }

        MeshBasedCellPopulation<3> cell_population(*p_mesh, cells);

    	// Make sure we have a Voronoi tessellation to begin with
        cell_population.CreateVoronoiTessellation();

        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputVoronoiData(true);
        cell_population.SetOutputCellMutationStates(true);
        cell_population.SetOutputCellProliferativeTypes(true);
        cell_population.SetOutputCellAncestors(true);

        // CryptSimulation3d(rCellPopulation, bool deleteCellPopulationAndForceCollection, bool initialiseCells, width, depth,
        // bool mRigidBoundaryConditions)
        CellBasedSimulation<3> simulator(cell_population);

		simulator.SetOutputDirectory("TestSimpleCube3d");
    	double normal_timestep = simulator.GetDt();		// This is 1/120 => every 30 seconds
    	double timestep_used = normal_timestep*0.5;
    	simulator.SetDt(timestep_used);
    	simulator.SetSamplingTimestepMultiple(240);		// Every hour
		simulator.SetEndTime(1.0);
        simulator.Solve();

        delete p_mesh;
    }

};

#endif /*TESTCELLBASEDSIMULATION3D_HPP_*/

