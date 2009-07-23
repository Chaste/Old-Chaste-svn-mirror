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
#ifndef TESTMESHBASEDTISSUE_HPP_
#define TESTMESHBASEDTISSUE_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "MeshBasedTissue.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModelCellsGenerator.hpp"
#include "AbstractCancerTestSuite.hpp"


class TestMeshBasedTissue : public AbstractCancerTestSuite
{
private:

    template<unsigned DIM>
    void TestSmallMeshBasedTissue(std::string meshFilename)
    {
        // Create a simple mesh
        TrianglesMeshReader<DIM,DIM> mesh_reader(meshFilename);
        MutableMesh<DIM,DIM> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node. Give each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<DIM> generator;
        generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create the tissue
        MeshBasedTissue<DIM> tissue(mesh, cells);

        TS_ASSERT_EQUALS(tissue.rGetMesh().GetNumNodes(), mesh.GetNumNodes());
        TS_ASSERT_EQUALS(tissue.rGetCells().size(), cells.size());

        unsigned counter = 0;
        for (typename AbstractTissue<DIM>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            // Test operator* and that cells are in sync
            TS_ASSERT_EQUALS(tissue.GetLocationIndexUsingCell(&(*cell_iter)), counter);

            // Test operator-> and that cells are in sync
            TS_ASSERT_DELTA(cell_iter->GetAge(), (double)counter, 1e-12);

            counter++;
        }

        TS_ASSERT_EQUALS(counter, tissue.GetNumRealCells());
    }


public:

    // Test construction, accessors and Iterator
    void TestSmallMeshBasedTissue1d2d3d() throw(Exception)
    {
        TestSmallMeshBasedTissue<1>("mesh/test/data/1D_0_to_1_10_elements");
        TestSmallMeshBasedTissue<2>("mesh/test/data/square_4_elements");
        TestSmallMeshBasedTissue<3>("mesh/test/data/cube_136_elements");
    }

    void TestValidateMeshBasedTissue()
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node apart from one.
        // Give each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<mesh.GetNumNodes()-1; i++)
        {
            AbstractCellCycleModel *p_cell_cycle_model = new FixedDurationGenerationBasedCellCycleModel();
            TissueCell cell(STEM, HEALTHY, p_cell_cycle_model);
            double birth_time = 0.0 - i;
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Fails as no cell corresponding to node 0
        TS_ASSERT_THROWS_ANYTHING(MeshBasedTissue<2> tissue2(mesh, cells));

        // Add another cell
        AbstractCellCycleModel *p_cell_cycle_model = new FixedDurationGenerationBasedCellCycleModel();
        TissueCell cell(STEM, HEALTHY, p_cell_cycle_model);
        double birth_time = -4.0;
        cell.SetBirthTime(birth_time);
        cells.push_back(cell);

        TS_ASSERT_THROWS_NOTHING(MeshBasedTissue<2> tissue2(mesh, cells));
    }

    void TestCreateCellPair()
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node. Give each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Give cells 0 and 1 specific mutations to enable later testing
        cells[0].SetMutationState(LABELLED);
        cells[1].SetMutationState(APC_ONE_HIT);

        // Create tissue
        MeshBasedTissue<2> tissue(mesh, cells);

        // Create cell pair
        std::set<TissueCell*> cell_pair = tissue.CreateCellPair(cells[0], cells[1]);

        // Check the cell pair was created correctly
        std::set<TissueCell*>::iterator cell_pair_iter = cell_pair.begin();

        TissueCell *p_cell0 = *cell_pair_iter;
        TS_ASSERT_EQUALS(p_cell0->GetMutationState(), LABELLED);

        ++cell_pair_iter;
        TissueCell *p_cell1 = *cell_pair_iter;
        TS_ASSERT_EQUALS(p_cell1->GetMutationState(), APC_ONE_HIT);
    }

    void TestAreaBasedDampingConstant()
    {
        // Create a simple mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0u, false);

        MutableMesh<2,2> *p_mesh = generator.GetMesh();

        // Set up cells, one for each node. Give each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());
        cells[9].SetMutationState(APC_TWO_HIT);

        MeshBasedTissue<2> tissue(*p_mesh, cells);

        TS_ASSERT_EQUALS(tissue.UseAreaBasedDampingConstant(), false);

        double damping_const = tissue.GetDampingConstant(8);

        TS_ASSERT_DELTA(damping_const, TissueConfig::Instance()->GetDampingConstantNormal(), 1e-6);

        double mutant_damping_const = tissue.GetDampingConstant(9);

        TS_ASSERT_DELTA(mutant_damping_const, TissueConfig::Instance()->GetDampingConstantMutant(), 1e-6);

        tissue.SetAreaBasedDampingConstant(true);

        TS_ASSERT_EQUALS(tissue.UseAreaBasedDampingConstant(), true);

        // Note that this method is usually called by TissueSimulation::Solve()
        tissue.CreateVoronoiTessellation();

        double area_based_damping_const = tissue.GetDampingConstant(8);

        // Since the tissue is in equilibrium, we should get the same damping constant as before
        TS_ASSERT_DELTA(area_based_damping_const, TissueConfig::Instance()->GetDampingConstantNormal(), 1e-6);
    }

    void TestSetNodeAndAddCell()
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a tissue, with no ghost nodes at the moment
        MeshBasedTissue<2> tissue(mesh, cells);

        //////////////////
        // Test set node
        //////////////////

        // Move node 0 by a small amount
        AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
        c_vector<double,2> new_location = tissue.GetLocationOfCellCentre(&(*cell_iter));
        new_location[0] += 1e-2;
        new_location[1] += 1e-2;
        ChastePoint<2> new_location_point(new_location);
        tissue.SetNode(tissue.GetLocationIndexUsingCell(&(*cell_iter)), new_location_point);

        TS_ASSERT_DELTA(mesh.GetNode(0)->rGetLocation()[0], new_location[0], 1e-12);
        TS_ASSERT_DELTA(mesh.GetNode(0)->rGetLocation()[1], new_location[1], 1e-12);

        //////////////////
        // Test add cell
        //////////////////
        unsigned old_num_nodes = mesh.GetNumNodes();
        unsigned old_num_cells = tissue.rGetCells().size();

        // Create a new cell, DON'T set the node index, set birth time=-1
        TissueCell cell(STEM, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
        cell.SetBirthTime(-1);
        c_vector<double,2> new_cell_location;
        new_cell_location[0] = 2;
        new_cell_location[1] = 2;

        tissue.AddCell(cell, new_cell_location);

        // Tissue should have updated mesh and cells
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), old_num_nodes+1);
        TS_ASSERT_EQUALS(tissue.rGetCells().size(), old_num_cells+1);
        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), old_num_nodes+1);

        // Same test via tissue class
        TS_ASSERT_EQUALS(tissue.rGetMesh().GetNumNodes(), old_num_nodes+1);
        TS_ASSERT_EQUALS(tissue.rGetCells().size(), old_num_cells+1);

        // Check the location of the new node
        TS_ASSERT_DELTA(mesh.GetNode(old_num_nodes)->rGetLocation()[0], 2.0, 1e-12);
        TS_ASSERT_DELTA(mesh.GetNode(old_num_nodes)->rGetLocation()[1], 2.0, 1e-12);

        // Check the index of the new cell
        TissueCell& new_cell = tissue.rGetCells().back();
        TS_ASSERT_EQUALS(tissue.GetLocationIndexUsingCell(&new_cell), old_num_nodes);
    }

    void TestRemoveDeadCellsAndUpdate()
    {
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells
        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
        cells[27].StartApoptosis();

        // Create a tissue without ghost nodes
        MeshBasedTissue<2> tissue(mesh,cells);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 81u);
        TS_ASSERT_EQUALS(tissue.rGetCells().size(), 81u);
        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), 81u);

        p_simulation_time->IncrementTimeOneStep();

        unsigned num_removed = tissue.RemoveDeadCells();

        TS_ASSERT_EQUALS(num_removed, 1u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 80u);
        TS_ASSERT_EQUALS(tissue.rGetCells().size(), 80u);
        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), 80u);
        TS_ASSERT_DIFFERS(tissue.rGetCells().size(), cells.size()); // Tissue now copies cells

        tissue.Update();

        // For coverage
        NodeMap map(mesh.GetNumAllNodes());
        map.ResetToIdentity();
        tissue.UpdateGhostNodesAfterReMesh(map);

        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), 80u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh.GetNumAllNodes());

        for (unsigned i=0; i<mesh.GetNumAllNodes(); i++)
        {
            TS_ASSERT_EQUALS(tissue.IsGhostNode(i), false);
        }

        // Finally, check the cells node indices have updated

        // We expect the cell node indices to be {0,11,...,79}
        std::set<unsigned> expected_node_indices;
        for (unsigned i=0; i<tissue.GetNumRealCells(); i++)
        {
            expected_node_indices.insert(i);
        }

        // Get actual cell node indices
        std::set<unsigned> node_indices;

        for (AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            // Record node index corresponding to cell
            unsigned node_index = tissue.GetLocationIndexUsingCell(&(*cell_iter));
            node_indices.insert(node_index);
        }

        TS_ASSERT_EQUALS(node_indices, expected_node_indices);
    }

    void TestUpdateNodeLocations()
    {
        // Test MeshBasedTissue::UpdateNodeLocations()

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a tissue, with no ghost nodes at the moment
        MeshBasedTissue<2> tissue(mesh, cells);

        // Make up some forces
        std::vector<c_vector<double, 2> > old_posns(tissue.GetNumNodes());
        std::vector<c_vector<double, 2> > forces_on_nodes(tissue.GetNumNodes());

        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
            old_posns[i][0] = tissue.GetNode(i)->rGetLocation()[0];
            old_posns[i][1] = tissue.GetNode(i)->rGetLocation()[1];

            forces_on_nodes[i][0] = i*0.01;
            forces_on_nodes[i][1] = 2*i*0.01;
        }

        // Call method
        double time_step = 0.01;
        tissue.UpdateNodeLocations(forces_on_nodes, time_step);

        // Check that node locations were correctly updated
        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA(tissue.GetNode(i)->rGetLocation()[0], old_posns[i][0] +   i*0.01*0.01, 1e-9);
            TS_ASSERT_DELTA(tissue.GetNode(i)->rGetLocation()[1], old_posns[i][1] + 2*i*0.01*0.01, 1e-9);
        }
    }

    void TestOutputWriters()
    {
        // Resetting the Maximum cell Id to zero (to account for previous tests)
        TissueCell::ResetMaxCellId();

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells
        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
        cells[0].SetCellType(APOPTOTIC); // coverage

        // Create tissue
        MeshBasedTissue<2> tissue(mesh, cells);

        // Test set methods

        TissueConfig::Instance()->SetOutputVoronoiData(true);
        TissueConfig::Instance()->SetOutputTissueAreas(true);
        TissueConfig::Instance()->SetOutputCellAreas(true);

        // This method is usually called by Update()
        tissue.CreateVoronoiTessellation();

        std::string output_directory = "TestTissueWriters";
        OutputFileHandler output_file_handler(output_directory, false);

		TissueConfig::Instance()->SetOutputCellMutationStates(true);
		TissueConfig::Instance()->SetOutputCellTypes(true);
		TissueConfig::Instance()->SetOutputCellAges(true);

        TS_ASSERT_THROWS_NOTHING(tissue.CreateOutputFiles(output_directory, false));

        tissue.WriteResultsToFiles();

        TS_ASSERT_THROWS_NOTHING(tissue.CloseOutputFiles());

        // Compare output with saved files of what they should look like
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizelements   cancer/test/data/TestTissueWriters/results.vizelements").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.viznodes      cancer/test/data/TestTissueWriters/results.viznodes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizcelltypes  cancer/test/data/TestTissueWriters/results.vizcelltypes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "tissueareas.dat       cancer/test/data/TestTissueWriters/tissueareas.dat").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellareas.dat       cancer/test/data/TestTissueWriters/cellareas.dat").c_str()), 0);

        // Test the GetCellMutationStateCount function: there should only be healthy cells
        c_vector<unsigned,5> cell_mutation_states = tissue.GetCellMutationStateCount();
        TS_ASSERT_EQUALS(cell_mutation_states[0], tissue.GetNumRealCells());
        for (unsigned i=1; i<NUM_CELL_MUTATION_STATES; i++)
        {
            TS_ASSERT_EQUALS(cell_mutation_states[i], 0u);
        }

        // Test the GetCellTypeCount function - we should have 4 stem cells and 1 dead cell (for coverage)
        c_vector<unsigned,5> cell_types = tissue.GetCellTypeCount();
        TS_ASSERT_EQUALS(cell_types[0], 4u);
        TS_ASSERT_EQUALS(cell_types[1], 0u);
        TS_ASSERT_EQUALS(cell_types[2], 0u);
        TS_ASSERT_EQUALS(cell_types[3], 1u);
    }

    void TestGetLocationOfCellCentreAndGetNodeCorrespondingToCell() throw (Exception)
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node. Give each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create the tissue
        MeshBasedTissue<2> tissue(mesh, cells);

        // Loop over nodes
        for (AbstractTissue<2>::Iterator cell_iter=tissue.Begin();
             cell_iter!=tissue.End();
             ++cell_iter)
        {
            // Record node location
            c_vector<double,2> node_location = tissue.GetLocationOfCellCentre(&(*cell_iter));

            // Test GetLocationOfCellCentre()
            TS_ASSERT_DELTA(node_location[0], tissue.GetLocationOfCellCentre(&(*cell_iter))[0], 1e-9);
            TS_ASSERT_DELTA(node_location[1], tissue.GetLocationOfCellCentre(&(*cell_iter))[1], 1e-9);
        }
    }

    // This test checks that the cells and nodes are correctly archived.
    void TestArchivingMeshBasedTissue() throw (Exception)
    {
        OutputFileHandler handler("archive",false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "mesh_based_tissue.arch";
        ArchiveLocationInfo::SetMeshPathname(handler.GetOutputDirectoryFullPath(),"mesh_based_tissue_mesh");

        std::vector<c_vector<double,2> > cell_locations;

        // Archive a tissue
        {
            // Need to set up time
            unsigned num_steps=10;
            SimulationTime *p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);

            // Create a simple mesh
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
            MutableMesh<2,2> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            // Set up cells, one for each node. Give each a birth time of -node_index,
            // so the age = node_index
            std::vector<TissueCell> cells;
            FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
            cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

            // Create the tissue
            MeshBasedTissue<2>* const p_tissue = new MeshBasedTissue<2>(mesh, cells);

            // Cells have been given birth times of 0, -1, -2, -3, -4.
            // loop over them to run to time 0.0;
            for (AbstractTissue<2>::Iterator cell_iter=p_tissue->Begin();
                 cell_iter!=p_tissue->End();
                 ++cell_iter)
            {
                cell_iter->ReadyToDivide();
                cell_locations.push_back(p_tissue->GetLocationOfCellCentre(&(*cell_iter)));
            }

            p_tissue->MarkSpring(p_tissue->rGetCellUsingLocationIndex(0), p_tissue->rGetCellUsingLocationIndex(1));

            // Set area-based viscosity
            p_tissue->SetAreaBasedDampingConstant(true);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Write the tissue to the archive
            output_arch << static_cast<const SimulationTime&> (*p_simulation_time);
            output_arch << p_tissue;
            SimulationTime::Destroy();
            delete p_tissue;
        }

        // Restore tissue
        {
            // Need to set up time
            unsigned num_steps=10;
            SimulationTime *p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);
            p_simulation_time->IncrementTimeOneStep();

            MeshBasedTissue<2> *p_tissue;

            // Restore the tissue
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> *p_simulation_time;

            input_arch >> p_tissue;

            // Cells have been given birth times of 0, -1, -2, -3, -4.
            // this checks that individual cells and their models are archived.
            unsigned counter = 0u;
            for (AbstractTissue<2>::Iterator cell_iter=p_tissue->Begin();
                 cell_iter!=p_tissue->End();
                 ++cell_iter)
            {
                TS_ASSERT_DELTA(cell_iter->GetAge(),(double)(counter),1e-7);
                TS_ASSERT_DELTA(p_tissue->GetLocationOfCellCentre(&(*cell_iter))[0],cell_locations[counter][0],1e-9);
                TS_ASSERT_DELTA(p_tissue->GetLocationOfCellCentre(&(*cell_iter))[1],cell_locations[counter][1],1e-9);
                counter++;
            }

            // Check the marked spring
            TS_ASSERT(p_tissue->IsMarkedSpring(p_tissue->rGetCellUsingLocationIndex(0), p_tissue->rGetCellUsingLocationIndex(1)));

            // Check the simulation time has been restored (through the cell)
            TS_ASSERT_EQUALS(p_simulation_time->GetTime(), 0.0);

            // Check the tissue has been restored
            TS_ASSERT_EQUALS(p_tissue->rGetCells().size(),5u);

            // Check area-based viscosity is still true
            TS_ASSERT_EQUALS(p_tissue->UseAreaBasedDampingConstant(), true);

            TS_ASSERT_EQUALS(p_tissue->rGetMesh().GetNumNodes(),5u);

            delete p_tissue;
        }
    }

    void TestSpringMarking()
    {
        // Create a small tissue
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0, 0.5));
        nodes.push_back(new Node<2>(1, false, 1, 0));
        nodes.push_back(new Node<2>(2, false, 1, 1));
        nodes.push_back(new Node<2>(3, false, 2, 0.5));
        nodes.push_back(new Node<2>(4, false, 2, 1.5));

        MutableMesh<2,2> mesh(nodes);

        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        MeshBasedTissue<2> tissue(mesh, cells);

        // Mark some springs
        tissue.MarkSpring(tissue.rGetCellUsingLocationIndex(1), tissue.rGetCellUsingLocationIndex(2));

        // Unmark and re-mark spring (for coverage)
        tissue.UnmarkSpring(tissue.rGetCellUsingLocationIndex(1), tissue.rGetCellUsingLocationIndex(2));
        tissue.MarkSpring(tissue.rGetCellUsingLocationIndex(1), tissue.rGetCellUsingLocationIndex(2));

        tissue.MarkSpring(tissue.rGetCellUsingLocationIndex(3), tissue.rGetCellUsingLocationIndex(4));

        // Check if springs are marked
        TS_ASSERT(tissue.IsMarkedSpring(tissue.rGetCellUsingLocationIndex(1), tissue.rGetCellUsingLocationIndex(2)));
        TS_ASSERT(tissue.IsMarkedSpring(tissue.rGetCellUsingLocationIndex(3), tissue.rGetCellUsingLocationIndex(4)));

        TS_ASSERT(!tissue.IsMarkedSpring(tissue.rGetCellUsingLocationIndex(1), tissue.rGetCellUsingLocationIndex(4)));
        TS_ASSERT(!tissue.IsMarkedSpring(tissue.rGetCellUsingLocationIndex(0), tissue.rGetCellUsingLocationIndex(2)));

        // Delete cell 4
        tissue.rGetCellUsingLocationIndex(4).Kill();
        tissue.RemoveDeadCells();

        // Check springs with non-deleted cells are still marked
        TS_ASSERT(tissue.IsMarkedSpring(tissue.rGetCellUsingLocationIndex(1), tissue.rGetCellUsingLocationIndex(2)));
        tissue.CheckTissueCellPointers();

        // Move cell 2
        AbstractTissue<2>::Iterator it = tissue.Begin();
        ++it;
        ++it;
        TS_ASSERT_EQUALS(tissue.GetLocationIndexUsingCell(&(*it)), 2u);
        ChastePoint<2> new_location(1, 10);
        tissue.SetNode(tissue.GetLocationIndexUsingCell(&(*it)), new_location);

        // Update tissue
        tissue.Update();

        tissue.CheckTissueCellPointers();

        // Check there is no marked spring between nodes 1 & 2
        TS_ASSERT(!tissue.IsMarkedSpring(tissue.rGetCellUsingLocationIndex(1), tissue.rGetCellUsingLocationIndex(2)));
    }

    void TestSettingCellAncestors() throw (Exception)
    {
        // Create a small tissue
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0, 0.5));
        nodes.push_back(new Node<2>(1, false, 1, 0));
        nodes.push_back(new Node<2>(2, false, 1, 1));
        nodes.push_back(new Node<2>(3, false, 2, 0.5));
        nodes.push_back(new Node<2>(4, false, 2, 1.5));

        MutableMesh<2,2> mesh(nodes);

        std::vector<TissueCell> cells;
        FixedDurationGenerationBasedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        MeshBasedTissue<2> tissue(mesh, cells);

        // Test that the tissue makes all cells fix the node index as ancestor
        tissue.SetCellAncestorsToNodeIndices();

        unsigned counter = 0;
        for (AbstractTissue<2>::Iterator cell_iter=tissue.Begin();
             cell_iter!=tissue.End();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS(cell_iter->GetAncestor(), tissue.GetLocationIndexUsingCell(&(*cell_iter)));
            counter ++;
        }
        TS_ASSERT_EQUALS(counter, 5u);

        // Test that we can recover remaining number of ancestors
        std::set<unsigned> remaining_ancestors = tissue.GetCellAncestors();
        TS_ASSERT_EQUALS(remaining_ancestors.size(), 5u);

        // Test that the set correctly represents a monoclonal population
        for (AbstractTissue<2>::Iterator cell_iter=tissue.Begin();
             cell_iter!=tissue.End();
             ++cell_iter)
        {
            // Set all cells to have the same ancestor
            cell_iter->SetAncestor(1u);
        }
        remaining_ancestors = tissue.GetCellAncestors();
        TS_ASSERT_EQUALS(remaining_ancestors.size(), 1u);
    }

};


#endif /*TESTMESHBASEDTISSUE_HPP_*/
