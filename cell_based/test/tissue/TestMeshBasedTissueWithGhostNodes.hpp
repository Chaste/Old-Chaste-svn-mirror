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
#ifndef TESTMESHBASEDTISSUEWITHGHOSTNODES_HPP_
#define TESTMESHBASEDTISSUEWITHGHOSTNODES_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "MeshBasedTissueWithGhostNodes.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CryptCellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"

class TestMeshBasedTissueWithGhostNodes : public AbstractCellBasedTestSuite
{
public:
    /*
     * Here we set up a test with 5 nodes, make a cell for each.
     * We then set cell 0 to be associated with node 1 instead of node 0
     * Validate throws an exception.
     * We then set node 0 to be a ghost node
     * Validate passes.
     */
    void TestValidateMeshBasedTissueWithGhostNodes()
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node apart from one.
        // Give each a birth time of -node_index, so the age = node_index
        std::vector<TissueCellPtr> cells;
        std::vector<unsigned> cell_location_indices;
        boost::shared_ptr<AbstractCellProperty> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<mesh.GetNumNodes()-1; i++)
        {
            AbstractCellCycleModel* p_cell_cycle_model = new FixedDurationGenerationBasedCellCycleModel();
            p_cell_cycle_model->SetCellProliferativeType(STEM);

            TissueCellPtr p_cell(new TissueCell(p_state, p_cell_cycle_model));
            double birth_time = 0.0 - i;
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
            cell_location_indices.push_back(i);
        }

        // Fails as the tissue constructor is not given the location indices
        // corresponding to real cells, so cannot work out which nodes are
        // ghost nodes
        std::vector<TissueCellPtr> cells_copy(cells);
        TS_ASSERT_THROWS_THIS(MeshBasedTissueWithGhostNodes<2> dodgy_tissue(mesh, cells_copy),
                "Node 4 does not appear to be a ghost node or have a cell associated with it");

        // Passes as the tissue constructor automatically works out which
        // cells are ghost nodes using the mesh and cell_location_indices
        MeshBasedTissueWithGhostNodes<2> tissue(mesh, cells, cell_location_indices);

        // Here we set the ghost nodes to what they already are
        std::set<unsigned> ghost_node_indices;
        ghost_node_indices.insert(mesh.GetNumNodes()-1u);
        tissue.SetGhostNodes(ghost_node_indices);

        // So validate passes at the moment
        tissue.Validate();

        // Test GetCellUsingLocationIndex()

        TS_ASSERT_THROWS_NOTHING(tissue.GetCellUsingLocationIndex(0)); // real cell
        TS_ASSERT_THROWS_THIS(tissue.GetCellUsingLocationIndex(mesh.GetNumNodes()-1u),"Location index input argument does not correspond to a TissueCell"); // ghost node

        // Now we label a real cell's node as a ghost
        ghost_node_indices.insert(1u);

        // Validate detects this inconsistency
        TS_ASSERT_THROWS_THIS(tissue.SetGhostNodes(ghost_node_indices),"Node 1 is labelled as a ghost node and has a cell attached");
    }

    // Test with ghost nodes, checking that the Iterator doesn't loop over ghost nodes
    void TestMeshBasedTissueWithGhostNodesSetup() throw(Exception)
    {
        unsigned num_cells_depth = 11;
        unsigned num_cells_width = 6;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2, false);

        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<TissueCellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel,2> cells_generator;
        cells_generator.GenerateGivenLocationIndices(cells, location_indices);

        // Create a tissue
        MeshBasedTissueWithGhostNodes<2> tissue(*p_mesh, cells, location_indices);

        // Create a set of node indices corresponding to ghost nodes
        std::set<unsigned> node_indices;
        std::set<unsigned> location_indices_set;
        std::set<unsigned> ghost_node_indices;

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            node_indices.insert(p_mesh->GetNode(i)->GetIndex());
        }
        for (unsigned i=0; i<location_indices.size(); i++)
        {
            location_indices_set.insert(location_indices[i]);
        }

        std::set_difference(node_indices.begin(), node_indices.end(),
                            location_indices_set.begin(), location_indices_set.end(),
                            std::inserter(ghost_node_indices, ghost_node_indices.begin()));

        std::vector<bool> is_ghost_node(p_mesh->GetNumNodes(), false);
        for (std::set<unsigned>::iterator it=ghost_node_indices.begin();
             it!=ghost_node_indices.end();
             it++)
        {
            is_ghost_node[*it] = true;
        }

        TS_ASSERT_EQUALS(tissue.rGetGhostNodes(), is_ghost_node);

        // Test the GetGhostNodeIndices method
        std::set<unsigned> ghost_node_indices2 = tissue.GetGhostNodeIndices();
        TS_ASSERT_EQUALS(ghost_node_indices, ghost_node_indices2);

        // Check the iterator doesn't loop over ghost nodes
        unsigned counter = 0;
        for (AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            unsigned node_index = tissue.GetLocationIndexUsingCell(*cell_iter);
            TS_ASSERT_EQUALS(is_ghost_node[node_index], false);
            counter++;
        }

        TS_ASSERT_EQUALS(counter, tissue.GetNumRealCells());

        // Check counter = num_nodes - num_ghost_nodes
        TS_ASSERT_EQUALS(counter + ghost_node_indices.size(), p_mesh->GetNumNodes());

        TS_ASSERT_EQUALS(tissue.rGetGhostNodes().size(), p_mesh->GetNumNodes());
    }


    void TestAreaBasedVisocityOnAPeriodicMesh() throw (Exception)
    {
        // Set up the simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        unsigned cells_across = 3;
        unsigned cells_up = 3;
        unsigned thickness_of_ghost_layer = 2;

        double scale_factor = 1.2;
        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, true, scale_factor);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<TissueCellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true); // true = mature cells

        MeshBasedTissueWithGhostNodes<2> tissue(*p_mesh, cells, location_indices);

        GeneralisedLinearSpringForce<2> linear_force;

        // It seems quite difficult to test this on a periodic mesh,
        // so just check the areas of all the cells are correct.

        tissue.CreateVoronoiTessellation();

        TS_ASSERT_EQUALS(tissue.GetVoronoiTessellation()->GetNumElements(), p_mesh->GetNumNodes());

        for (AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            unsigned node_index = tissue.GetLocationIndexUsingCell(*cell_iter);
            double area = tissue.GetVolumeOfVoronoiElement(node_index);
            TS_ASSERT_DELTA(area, sqrt(3)*scale_factor*scale_factor/2, 1e-6);
        }
    }


    void TestRemoveDeadCellsAndReMeshWithGhostNodes()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create vector of cell location indices
        std::vector<unsigned> cell_location_indices;
        for (unsigned i=10; i<mesh.GetNumNodes(); i++)
        {
            if (i!=80)
            {
                cell_location_indices.push_back(i);
            }
        }

        // Set up cells
        std::vector<TissueCellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, cell_location_indices.size());
        cells[27]->StartApoptosis();

        // Create a tissue, with some random ghost nodes
        MeshBasedTissueWithGhostNodes<2> tissue_with_ghost_nodes(mesh, cells, cell_location_indices);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 81u);

        // Num real cells should be num_nodes (81) - num_ghosts (11) = 70
        TS_ASSERT_EQUALS(tissue_with_ghost_nodes.GetNumRealCells(), 70u);

        p_simulation_time->IncrementTimeOneStep();

        unsigned num_removed_with_ghost_nodes = tissue_with_ghost_nodes.RemoveDeadCells();

        TS_ASSERT_EQUALS(num_removed_with_ghost_nodes, 1u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 80u);
        TS_ASSERT_DIFFERS(tissue_with_ghost_nodes.rGetCells().size(), cells.size()); // Tissue now copies cells

        // Num real cells should be num_nodes (81) - num_ghosts (11) - 1 deleted node = 69
        TS_ASSERT_EQUALS(tissue_with_ghost_nodes.GetNumRealCells(), 69u);
        TS_ASSERT_EQUALS(tissue_with_ghost_nodes.rGetGhostNodes().size(), mesh.GetNumAllNodes());

        tissue_with_ghost_nodes.Update();

        // For coverage
        NodeMap map(mesh.GetNumAllNodes());
        map.ResetToIdentity();
        tissue_with_ghost_nodes.UpdateGhostNodesAfterReMesh(map);

        // Num real cells should be new_num_nodes (80) - num_ghosts (11)
        TS_ASSERT_EQUALS(tissue_with_ghost_nodes.GetNumRealCells(), 69u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh.GetNumAllNodes());
        TS_ASSERT_EQUALS(tissue_with_ghost_nodes.rGetGhostNodes().size(), mesh.GetNumNodes());

        // Nodes 0-9 should not been renumbered so are still ghost nodes.
        // the ghost node at node 80 is now at 79 as node 27 was deleted..
        for (unsigned i=0; i<mesh.GetNumAllNodes(); i++)
        {
            // True (ie should be a ghost) if i<10 or i==79, else false
            TS_ASSERT_EQUALS(tissue_with_ghost_nodes.IsGhostNode(i), ((i<10)||(i==79)));
        }

        // Finally, check the cells node indices have updated

        // We expect the cell node indices to be {10,11,...,79}
        std::set<unsigned> expected_node_indices;
        for (unsigned i=0; i<tissue_with_ghost_nodes.GetNumRealCells(); i++)
        {
            expected_node_indices.insert(i+10);
        }

        // Get actual cell node indices
        std::set<unsigned> node_indices_with_ghost_nodes;

        for (AbstractTissue<2>::Iterator cell_iter = tissue_with_ghost_nodes.Begin();
             cell_iter != tissue_with_ghost_nodes.End();
             ++cell_iter)
        {
            // Record node index corresponding to cell
            unsigned node_index_with_ghost_nodes = tissue_with_ghost_nodes.GetLocationIndexUsingCell(*cell_iter);
            node_indices_with_ghost_nodes.insert(node_index_with_ghost_nodes);
        }

        TS_ASSERT_EQUALS(node_indices_with_ghost_nodes, expected_node_indices);
    }

    void TestAddAndRemoveAndAddWithOutUpdate()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create vector of cell location indices
        std::vector<unsigned> cell_location_indices;
        for (unsigned i=10; i<mesh.GetNumNodes(); i++)
        {
            if (i!=80)
            {
                cell_location_indices.push_back(i);
            }
        }

        // Set up cells
        std::vector<TissueCellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, cell_location_indices.size());
        cells[27]->StartApoptosis();

        // Create a tissue, with some random ghost nodes
        MeshBasedTissueWithGhostNodes<2> tissue(mesh, cells, cell_location_indices);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 81u);
        TS_ASSERT_EQUALS(tissue.rGetCells().size(), 70u);

        boost::shared_ptr<AbstractCellProperty> p_state(new WildTypeCellMutationState);

        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
        p_model->SetCellProliferativeType(STEM);
        TissueCellPtr p_new_cell(new TissueCell(p_state, p_model));
        p_new_cell->SetBirthTime(0);

        c_vector<double,2> new_location;
        new_location[0] = 0.3433453454443;
        new_location[0] = 0.3435346344234;
        tissue.AddCell(p_new_cell, new_location);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 82u);
        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), 71u);

        p_simulation_time->IncrementTimeOneStep();

        unsigned num_removed = tissue.RemoveDeadCells();
        TS_ASSERT_EQUALS(num_removed, 1u);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 81u);
        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), 70u);

        FixedDurationGenerationBasedCellCycleModel* p_model2 = new FixedDurationGenerationBasedCellCycleModel();
        p_model2->SetCellProliferativeType(STEM);
        TissueCellPtr p_new_cell2(new TissueCell(p_state, p_model2));
        p_new_cell2->SetBirthTime(0);

        c_vector<double,2> new_location2;
        new_location2[0] = 0.6433453454443;
        new_location2[0] = 0.6435346344234;
        tissue.AddCell(p_new_cell2, new_location2);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 82u);
        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), 71u);
    }


    void TestUpdateNodeLocations() throw(Exception)
    {
        HoneycombMeshGenerator generator(3, 3, 1, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        std::vector<TissueCellPtr> cells2;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator2;
        cells_generator2.Generate(cells2, p_mesh, location_indices, true);

        MeshBasedTissueWithGhostNodes<2> tissue2(*p_mesh, cells2, location_indices);

        // Make up some forces
        std::vector<c_vector<double, 2> > old_posns2(tissue2.GetNumNodes());
        std::vector<c_vector<double, 2> > forces_on_nodes2(tissue2.GetNumNodes());

        for (unsigned i=0; i<tissue2.GetNumNodes(); i++)
        {
            old_posns2[i][0] = tissue2.GetNode(i)->rGetLocation()[0];
            old_posns2[i][1] = tissue2.GetNode(i)->rGetLocation()[1];

            forces_on_nodes2[i][0] = i*0.01;
            forces_on_nodes2[i][1] = 2*i*0.01;
        }

        // Call method
        double time_step = 0.01;
        tissue2.UpdateNodeLocations(forces_on_nodes2, time_step);

        // Check that node locations were correctly updated
        for (AbstractTissue<2>::Iterator cell_iter = tissue2.Begin();
             cell_iter != tissue2.End();
             ++cell_iter)
        {
            unsigned i = tissue2.GetLocationIndexUsingCell(*cell_iter);
            TS_ASSERT_DELTA(tissue2.GetNode(i)->rGetLocation()[0], old_posns2[i][0] +   i*0.01*0.01, 1e-9);
            TS_ASSERT_DELTA(tissue2.GetNode(i)->rGetLocation()[1], old_posns2[i][1] + 2*i*0.01*0.01, 1e-9);
        }
    }



    void TestSpringIterator2d() throw(Exception)
    {
        // Set up expected results for the honeycombmesh created below
        // the following are the edges which do not contain a ghost node
        std::set < std::set < unsigned > > expected_node_pairs;
        unsigned expected_node_pairs_array[] = {5,6,
                                                5,9,
                                                5,10,
                                                6,10,
                                                9,10 };

        for (unsigned i=0; i<10; i=i+2)
        {
            std::set < unsigned > node_pair;
            node_pair.insert(expected_node_pairs_array[i]);
            node_pair.insert(expected_node_pairs_array[i+1]);
            expected_node_pairs.insert(node_pair);
        }

        // set up simple tissue with honeycomb mesh

        unsigned num_cells_depth = 2;
        unsigned num_cells_width = 2;
        unsigned thickness_of_ghosts = 1;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, thickness_of_ghosts, false);

        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<TissueCellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel,2> cells_generator;
        cells_generator.GenerateGivenLocationIndices(cells, location_indices);

        // Create a tissue
        MeshBasedTissueWithGhostNodes<2> tissue(*p_mesh, cells, location_indices);

        // Check that we can iterate over the set of springs
        std::set< std::set< unsigned > > springs_visited;

        for (MeshBasedTissue<2>::SpringIterator spring_iterator=tissue.SpringsBegin();
             spring_iterator!=tissue.SpringsEnd();
             ++spring_iterator)
        {
            std::set<unsigned> node_pair;
            node_pair.insert(spring_iterator.GetNodeA()->GetIndex());
            node_pair.insert(spring_iterator.GetNodeB()->GetIndex());

            TS_ASSERT_EQUALS(springs_visited.find(node_pair), springs_visited.end());
            springs_visited.insert(node_pair);

            TS_ASSERT_EQUALS(tissue.GetLocationIndexUsingCell(spring_iterator.GetCellA()),
                             spring_iterator.GetNodeA()->GetIndex());

            TS_ASSERT_EQUALS(tissue.GetLocationIndexUsingCell(spring_iterator.GetCellB()),
                             spring_iterator.GetNodeB()->GetIndex());
        }

         TS_ASSERT_EQUALS(springs_visited, expected_node_pairs);
    }

    // 3d test with some ghost nodes
    void TestSpringIterator3d() throw(Exception)
    {
        // Create a simple mesh
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements");
        MutableMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create vector of cell location indices
        std::vector<unsigned> cell_location_indices;
        for (unsigned i=10; i<mesh.GetNumNodes(); i++)
        {
            cell_location_indices.push_back(i);
        }

        // Set up cells
        std::vector<TissueCellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> generator;
        generator.GenerateBasic(cells, cell_location_indices.size());

        // Create a tissue, with no ghost nodes at the moment
        MeshBasedTissueWithGhostNodes<3> tissue(mesh, cells, cell_location_indices);

        // Check that we can iterate over the set of springs
        std::set< std::set< unsigned > > springs_visited;

        for (MeshBasedTissue<3>::SpringIterator spring_iterator=tissue.SpringsBegin();
             spring_iterator!=tissue.SpringsEnd();
             ++spring_iterator)
        {
            std::set<unsigned> node_pair;
            node_pair.insert(spring_iterator.GetNodeA()->GetIndex());
            node_pair.insert(spring_iterator.GetNodeB()->GetIndex());

            TS_ASSERT_EQUALS(springs_visited.find(node_pair), springs_visited.end());
            springs_visited.insert(node_pair);

            TS_ASSERT_EQUALS(tissue.GetLocationIndexUsingCell(spring_iterator.GetCellA()),
                             spring_iterator.GetNodeA()->GetIndex());

            TS_ASSERT_EQUALS(tissue.GetLocationIndexUsingCell(spring_iterator.GetCellB()),
                             spring_iterator.GetNodeB()->GetIndex());
        }

        // Set up expected node pairs
        std::set< std::set<unsigned> > expected_node_pairs;
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            Element<3,3>* p_element = mesh.GetElement(i);
            for (unsigned j=0; j<4; j++)
            {
                for (unsigned k=0; k<4; k++)
                {
                    unsigned node_A = p_element->GetNodeGlobalIndex(j);
                    unsigned node_B = p_element->GetNodeGlobalIndex(k);

                    // If nodeA or node_B are <10 they will have been labelled a ghost node above
                    if (node_A != node_B && node_A>=10 && node_B>=10)
                    {
                        std::set<unsigned> node_pair;
                        node_pair.insert(node_A);
                        node_pair.insert(node_B);

                        expected_node_pairs.insert(node_pair);
                    }
                }
            }
        }

        TS_ASSERT_EQUALS(springs_visited, expected_node_pairs);
    }

    void TestTissueWritersIn3dWithGhostNodes()
    {
        // Set up SimulationTime (needed if VTK is used)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Resetting the Maximum cell Id to zero (to account for previous tests)
        TissueCell::ResetMaxCellId();

        // Create a simple 3D mesh with some ghost nodes
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0,  true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1,  true,  1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(2,  true,  1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(3,  true,  0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(4,  false, 0.5, 0.5, 0.5));
        nodes.push_back(new Node<3>(5,  false, -1.0, -1.0, -1.0));
        nodes.push_back(new Node<3>(6,  false,  2.0, -1.0, -1.0));
        nodes.push_back(new Node<3>(7,  false,  2.0,  2.0, -1.0));
        nodes.push_back(new Node<3>(8,  false, -1.0,  2.0, -1.0));
        nodes.push_back(new Node<3>(9,  false, -1.0, -1.0,  2.0));
        nodes.push_back(new Node<3>(10, false,  2.0, -1.0,  2.0));
        nodes.push_back(new Node<3>(11, false,  2.0,  2.0,  2.0));
        nodes.push_back(new Node<3>(12, false, -1.0,  2.0,  2.0));
        MutableMesh<3,3> mesh(nodes);

        std::vector<unsigned> location_indices;
        for (unsigned index=0; index<5; index++)
        {
            location_indices.push_back(index);
        }

        // Set up cells
        std::vector<TissueCellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateGivenLocationIndices(cells, location_indices);

        cells[4]->AddCellProperty(CellPropertyRegistry::Instance()->Get<ApoptoticCellProperty>()); // coverage

        TS_ASSERT_EQUALS(cells[4]->HasCellProperty<ApoptoticCellProperty>(), true);

        // Create tissue
        MeshBasedTissueWithGhostNodes<3> tissue(mesh, cells, location_indices);

        TS_ASSERT_EQUALS(tissue.GetIdentifier(), "MeshBasedTissueWithGhostNodes<3>");

        // Test set methods
        tissue.SetOutputVoronoiData(true);
        tissue.SetOutputTissueVolumes(true);
        tissue.SetOutputCellVolumes(true);
        tissue.SetOutputCellAncestors(true);
        tissue.SetOutputCellMutationStates(true);
        tissue.SetOutputCellProliferativeTypes(true);
        tissue.SetOutputCellAges(true);
        tissue.SetOutputCellCyclePhases(true);

        // This method is usually called by Update()
        tissue.CreateVoronoiTessellation();

        std::string output_directory = "TestTissueWritersIn3dWithGhostNodes";
        OutputFileHandler output_file_handler(output_directory, false);

        tissue.CreateOutputFiles(output_directory, false);
        tissue.WriteResultsToFiles();
        tissue.CloseOutputFiles();

        // Compare output with saved files of what they should look like
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizelements   cell_based/test/data/TestTissueWritersIn3dWithGhostNodes/results.vizelements").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.viznodes      cell_based/test/data/TestTissueWritersIn3dWithGhostNodes/results.viznodes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizcelltypes  cell_based/test/data/TestTissueWritersIn3dWithGhostNodes/results.vizcelltypes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "tissueareas.dat       cell_based/test/data/TestTissueWritersIn3dWithGhostNodes/tissueareas.dat").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellareas.dat         cell_based/test/data/TestTissueWritersIn3dWithGhostNodes/cellareas.dat").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "voronoi.dat           cell_based/test/data/TestTissueWritersIn3dWithGhostNodes/voronoi.dat").c_str()), 0);

        // Test the GetCellMutationStateCount function: there should only be healthy cells
        std::vector<unsigned> cell_mutation_states = tissue.GetCellMutationStateCount();
        TS_ASSERT_EQUALS(cell_mutation_states.size(), 4u);
        TS_ASSERT_EQUALS(cell_mutation_states[0], 5u);
        TS_ASSERT_EQUALS(cell_mutation_states[1], 0u);
        TS_ASSERT_EQUALS(cell_mutation_states[2], 0u);
        TS_ASSERT_EQUALS(cell_mutation_states[3], 0u);

        // Test the GetCellProliferativeTypeCount function - we should have 4 stem cells and 1 dead cell (for coverage)
        std::vector<unsigned> cell_types = tissue.rGetCellProliferativeTypeCount();
        TS_ASSERT_EQUALS(cell_types.size(), 3u);
        TS_ASSERT_EQUALS(cell_types[0], 5u);
        TS_ASSERT_EQUALS(cell_types[1], 0u);
        TS_ASSERT_EQUALS(cell_types[2], 0u);

        //Test that the tissue parameters are output correctly
		out_stream ParameterFile = output_file_handler.OpenOutputFile("results.parameters");
		// Write tissue parameters to file
		tissue.OutputTissueParameters(ParameterFile);
		ParameterFile->close();

		// Compare output with saved files of what they should look like
		TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.parameters     	cell_based/test/data/TestTissueWritersIn3dWithGhostNodes/results.parameters").c_str()), 0);

    }

    void TestVoronoiAreasAndPerimetersWithGhostNodes() throw (Exception)
    {
        // Create a small honeycomb mesh surrounded by a single layer of ghost nodes
        HoneycombMeshGenerator generator(2, 2, 1, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create some cells
        std::vector<TissueCellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel,2> cells_generator;
        cells_generator.GenerateGivenLocationIndices(cells, location_indices);

        // Create a tissue
        MeshBasedTissueWithGhostNodes<2> tissue(*p_mesh, cells, location_indices);

        // Create Voronoi tessellation (normally done in a simulation)
        tissue.CreateVoronoiTessellation();

        // The Voronoi element corresponding to each real cell should be a regular hexagon
        for (AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            unsigned node_index = tissue.GetLocationIndexUsingCell(*cell_iter);
            if (!tissue.IsGhostNode(node_index))
            {
                TS_ASSERT_DELTA(tissue.GetVolumeOfVoronoiElement(node_index), sqrt(3)/2, 1e-4);
                TS_ASSERT_DELTA(tissue.GetSurfaceAreaOfVoronoiElement(node_index), 6/sqrt(3), 1e-4);
            }
        }
    }

    void TestVoronoiGhostNodeLabelling2d() throw (Exception)
    {
        // Set up the simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        // Create 2D mesh with ghost nodes
        MutableMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(6, 6);

        c_vector<double, 2> mesh_centre = zero_vector<double>(2);
        for (unsigned node_index=0; node_index<mesh.GetNumNodes(); node_index++)
        {
            mesh_centre += mesh.GetNode(node_index)->rGetLocation() / mesh.GetNumNodes();
        }

        // Set up cells by iterating through the nodes
        std::vector<TissueCellPtr> cells;
        std::vector<unsigned> location_indices;
        boost::shared_ptr<AbstractCellProperty> p_state(new WildTypeCellMutationState);

        // Loop over nodes
        for (unsigned node_index=0; node_index<mesh.GetNumNodes(); node_index++)
        {
            // If this node is sufficiently close to the centre of the mesh, then create a cell for it
            c_vector<double, 2> node_location = mesh.GetNode(node_index)->rGetLocation();
            if (node_location(0) <= 3)
            {
                FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
                p_model->SetCellProliferativeType(TRANSIT);

                TissueCellPtr p_cell(new TissueCell(p_state, p_model));
                p_cell->SetBirthTime(-1.0);

                cells.push_back(p_cell);
                location_indices.push_back(node_index);
            }
        }

        // Create tissue
        MeshBasedTissueWithGhostNodes<2> tissue(mesh, cells, location_indices);
        tissue.SetOutputVoronoiData(true);

        // Create Voronoi tessellation
        tissue.CreateVoronoiTessellation();

        // Check the correspondence between ghost nodes is correct
        for (VertexMesh<2,2>::VertexElementIterator elem_iter = tissue.GetVoronoiTessellation()->GetElementIteratorBegin();
             elem_iter != tissue.GetVoronoiTessellation()->GetElementIteratorEnd();
             ++elem_iter)
        {
            unsigned elem_index = elem_iter->GetIndex();
            unsigned node_index = tissue.GetVoronoiTessellation()->GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex(elem_index);

            c_vector<double, 2> node_location = tissue.GetNode(node_index)->rGetLocation();
            bool should_be_ghost_node = (node_location(0) > 3);

            TS_ASSERT_EQUALS(tissue.IsGhostNode(node_index), should_be_ghost_node);
        }
    }

    void TestVoronoiGhostNodeLabelling3d() throw (Exception)
    {
        // Set up the simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        // Create 3D mesh with ghost nodes
        MutableMesh<3,3> mesh;
        mesh.ConstructCuboid(6, 6, 6);

        c_vector<double, 3> mesh_centre = zero_vector<double>(3);
        for (unsigned node_index=0; node_index<mesh.GetNumNodes(); node_index++)
        {
            mesh_centre += mesh.GetNode(node_index)->rGetLocation() / mesh.GetNumNodes();
        }

        // Set up cells by iterating through the nodes
        std::vector<TissueCellPtr> cells;
        std::vector<unsigned> location_indices;
        boost::shared_ptr<AbstractCellProperty> p_state(new WildTypeCellMutationState);

        // Loop over nodes
        for (unsigned node_index=0; node_index<mesh.GetNumNodes(); node_index++)
        {
            // If this node is sufficiently close to the centre of the mesh, then create a cell for it
            c_vector<double, 3> node_location = mesh.GetNode(node_index)->rGetLocation();
            if (node_location(0) <= 3)
            {
                FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
                p_model->SetCellProliferativeType(TRANSIT);

                TissueCellPtr p_cell(new TissueCell(p_state, p_model));
                p_cell->SetBirthTime(-1.0);

                cells.push_back(p_cell);
                location_indices.push_back(node_index);
            }
        }

        // Create tissue
        MeshBasedTissueWithGhostNodes<3> tissue(mesh, cells, location_indices);
        tissue.SetOutputVoronoiData(true);

        // Create Voronoi tessellation
        tissue.CreateVoronoiTessellation();

        // Check the correspondence between ghost nodes is correct
        for (VertexMesh<3,3>::VertexElementIterator elem_iter = tissue.GetVoronoiTessellation()->GetElementIteratorBegin();
             elem_iter != tissue.GetVoronoiTessellation()->GetElementIteratorEnd();
             ++elem_iter)
        {
            unsigned elem_index = elem_iter->GetIndex();
            unsigned node_index = tissue.GetVoronoiTessellation()->GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex(elem_index);

            c_vector<double, 3> node_location = tissue.GetNode(node_index)->rGetLocation();
            bool should_be_ghost_node = (node_location(0) > 3);

            TS_ASSERT_EQUALS(tissue.IsGhostNode(node_index), should_be_ghost_node);
        }
    }
};


#endif /*TESTMESHBASEDTISSUEWITHGHOSTNODES_HPP_*/
