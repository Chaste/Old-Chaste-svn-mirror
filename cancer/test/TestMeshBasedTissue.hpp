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

#include "MeshBasedTissueWithGhostNodes.hpp"
#include "MeinekeInteractionForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "FixedCellCycleModelCellsGenerator.hpp"
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

        // Set up cells, one for each node. Get each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<DIM> generator;
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
        // Get each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<mesh.GetNumNodes()-1; i++)
        {
            AbstractCellCycleModel* p_cell_cycle_model = new FixedCellCycleModel();
            TissueCell cell(STEM, HEALTHY, p_cell_cycle_model);
            double birth_time = 0.0 - i;
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Fails as no cell corresponding to node 0
        TS_ASSERT_THROWS_ANYTHING(MeshBasedTissue<2> tissue2(mesh, cells));

        // Add another cell
        AbstractCellCycleModel* p_cell_cycle_model = new FixedCellCycleModel();
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

        // Set up cells, one for each node. Get each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
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

        TissueCell* p_cell0 = *cell_pair_iter;
        TS_ASSERT_EQUALS(p_cell0->GetMutationState(), LABELLED);

        ++cell_pair_iter;
        TissueCell* p_cell1 = *cell_pair_iter;
        TS_ASSERT_EQUALS(p_cell1->GetMutationState(), APC_ONE_HIT);
    }

    void TestAreaBasedDampingConstant()
    {
        // Create a simple mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0u, false);

        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells, one for each node. Get each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());
        cells[9].SetMutationState(APC_TWO_HIT);

        MeshBasedTissue<2> tissue(*p_mesh, cells);

        TS_ASSERT_EQUALS(tissue.UseAreaBasedDampingConstant(), false);

        double damping_const = tissue.GetDampingConstant(8);

        TS_ASSERT_DELTA(damping_const, CancerParameters::Instance()->GetDampingConstantNormal(), 1e-6);

        double mutant_damping_const = tissue.GetDampingConstant(9);

        TS_ASSERT_DELTA(mutant_damping_const, CancerParameters::Instance()->GetDampingConstantMutant(), 1e-6);

        tissue.SetAreaBasedDampingConstant(true);

        TS_ASSERT_EQUALS(tissue.UseAreaBasedDampingConstant(), true);

        // Note that this method is usually called by TissueSimulation::Solve()
        tissue.CreateVoronoiTessellation();

        double area_based_damping_const = tissue.GetDampingConstant(8);

        // Since the tissue is in equilibrium, we should get the same damping constant as before
        TS_ASSERT_DELTA(area_based_damping_const, CancerParameters::Instance()->GetDampingConstantNormal(), 1e-6);
    }

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
        // Get each a birth time of -node_index, so the age = node_index
        std::vector<TissueCell> cells;
        std::vector<unsigned> cell_location_indices;
        for (unsigned i=0; i<mesh.GetNumNodes()-1; i++)
        {
            AbstractCellCycleModel* p_cell_cycle_model = new FixedCellCycleModel();
            TissueCell cell(STEM, HEALTHY, p_cell_cycle_model);
            double birth_time = 0.0 - i;
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
            
            cell_location_indices.push_back(i);
        }

        // Fails as the tissue constructor is not given the location indices
        // corresponding to real cells, so cannot work out which nodes are
        // ghost nodes
        TS_ASSERT_THROWS_ANYTHING(MeshBasedTissueWithGhostNodes<2> tissue(mesh, cells));
        
        // Passes as the tissue constructor automatically works out which 
        // cells are ghost nodes using the mesh and cell_location_indices
        TS_ASSERT_THROWS_NOTHING(MeshBasedTissueWithGhostNodes<2> tissue(mesh, cells, cell_location_indices));
    }

    // Test with ghost nodes, checking that the Iterator doesn't loop over ghost nodes
    void TestMeshBasedTissueWithGhostNodes() throw(Exception)
    {
        unsigned num_cells_depth = 11;
        unsigned num_cells_width = 6;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2u, false);

        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells 
        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
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
            unsigned node_index = tissue.GetLocationIndexUsingCell(&(*cell_iter));
            TS_ASSERT_EQUALS(is_ghost_node[node_index], false);
            counter++;
        }

        TS_ASSERT_EQUALS(counter, tissue.GetNumRealCells());

        // Check counter = num_nodes - num_ghost_nodes
        TS_ASSERT_EQUALS(counter + ghost_node_indices.size(), p_mesh->GetNumNodes());

        TS_ASSERT_EQUALS(tissue.rGetGhostNodes().size(), p_mesh->GetNumNodes());
    }

    void TestSetNodeAndAddCell()
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a tissue, with no ghost nodes at the moment
        MeshBasedTissue<2> tissue(mesh, cells);

        //////////////////
        // Test set node
        //////////////////

        // Move node 0 by a small amount
        AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
        c_vector<double,2> new_location = tissue.GetLocationOfCell(&(*cell_iter));
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
        TissueCell cell(STEM, HEALTHY, new FixedCellCycleModel());
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

    void TestAreaBasedVisocityOnAPeriodicMesh() throw (Exception)
    {
        // Set up the simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

        unsigned cells_across = 3;
        unsigned cells_up = 3;
        unsigned thickness_of_ghost_layer = 2;

        double scale_factor = 1.2;
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, true, scale_factor);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateForCrypt(cells, *p_mesh, location_indices, true); // true = mature cells

        MeshBasedTissueWithGhostNodes<2> tissue(*p_mesh, cells, location_indices);

        MeinekeInteractionForce<2> meineke_force;

        // It seems quite difficult to test this on a periodic mesh,
        // so just check the areas of all the cells are correct.

        tissue.CreateVoronoiTessellation();

        for (AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            unsigned node_index = tissue.GetLocationIndexUsingCell(&(*cell_iter));
            VoronoiTessellation<2>& tess = tissue.rGetVoronoiTessellation();
            double area = tess.GetFaceArea(node_index);
            TS_ASSERT_DELTA(area, sqrt(3)*scale_factor*scale_factor/2, 1e-6);
        }
    }

    void TestRemoveDeadCellsAndUpdate()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells
        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
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
        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateBasic(cells, cell_location_indices.size());
        cells[27].StartApoptosis();

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
            unsigned node_index_with_ghost_nodes = tissue_with_ghost_nodes.GetLocationIndexUsingCell(&(*cell_iter));
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
        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateBasic(cells, cell_location_indices.size());
        cells[27].StartApoptosis();

        // Create a tissue, with some random ghost nodes
        MeshBasedTissueWithGhostNodes<2> tissue(mesh, cells, cell_location_indices);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 81u);
        TS_ASSERT_EQUALS(tissue.rGetCells().size(), 70u);

        TissueCell new_cell(STEM, HEALTHY, new FixedCellCycleModel());
        new_cell.SetBirthTime(0);

        c_vector<double,2> new_location;
        new_location[0] = 0.3433453454443;
        new_location[0] = 0.3435346344234;
        tissue.AddCell(new_cell, new_location);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 82u);
        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), 71u);

        p_simulation_time->IncrementTimeOneStep();

        unsigned num_removed = tissue.RemoveDeadCells();
        TS_ASSERT_EQUALS(num_removed, 1u);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 81u);
        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), 70u);

        TissueCell new_cell2(STEM, HEALTHY, new FixedCellCycleModel());
        new_cell2.SetBirthTime(0);

        c_vector<double,2> new_location2;
        new_location2[0] = 0.6433453454443;
        new_location2[0] = 0.6435346344234;
        tissue.AddCell(new_cell2, new_location2);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 82u);
        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), 71u);
    }

    void TestUpdateNodeLocations()
    {
        // Test MeshBasedTissue::UpdateNodeLocations()

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
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

        // Test MeshBasedTissueWithGhostNodes::UpdateNodeLocations()

        HoneycombMeshGenerator generator(3, 3, 1, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        std::vector<TissueCell> cells2;
        FixedCellCycleModelCellsGenerator<2> cells_generator2;
        cells_generator2.GenerateForCrypt(cells2, *p_mesh, location_indices, true);

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
        tissue2.UpdateNodeLocations(forces_on_nodes2, time_step);

        // Check that node locations were correctly updated
        for (AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            unsigned i = tissue2.GetLocationIndexUsingCell(&(*cell_iter));
            TS_ASSERT_DELTA(tissue2.GetNode(i)->rGetLocation()[0], old_posns2[i][0] +   i*0.01*0.01, 1e-9);
            TS_ASSERT_DELTA(tissue2.GetNode(i)->rGetLocation()[1], old_posns2[i][1] + 2*i*0.01*0.01, 1e-9);
        }
    }

    void TestOutputWriters()
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells
        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
        cells[0].SetCellType(APOPTOTIC); // coverage
        MeshBasedTissue<2> tissue(mesh,cells);
        
        // Create tissue
        tissue.SetWriteTissueAreas(true); // coverage

        std::string output_directory = "TestTissueWriters";
        OutputFileHandler output_file_handler(output_directory, false);

        TS_ASSERT_THROWS_NOTHING(tissue.CreateOutputFiles(output_directory, false, true, true, false, false, false));

        tissue.WriteResultsToFiles(true, true, false, false, false);

        TS_ASSERT_THROWS_NOTHING(tissue.CloseOutputFiles(true, true, false, false, false));

        // Compare output with saved files of what they should look like
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizelements   cancer/test/data/TestTissueWriters/results.vizelements").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.viznodes      cancer/test/data/TestTissueWriters/results.viznodes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizcelltypes  cancer/test/data/TestTissueWriters/results.vizcelltypes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "tissueareas.dat       cancer/test/data/TestTissueWriters/tissueareas.dat").c_str()), 0);

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
        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
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

            TS_ASSERT_EQUALS(tissue.GetLocationIndexUsingCell(&(spring_iterator.rGetCellA())),
                             spring_iterator.GetNodeA()->GetIndex());

            TS_ASSERT_EQUALS(tissue.GetLocationIndexUsingCell(&(spring_iterator.rGetCellB())),
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
        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<3> generator;
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

            TS_ASSERT_EQUALS(tissue.GetLocationIndexUsingCell(&(spring_iterator.rGetCellA())),
                             spring_iterator.GetNodeA()->GetIndex());

            TS_ASSERT_EQUALS(tissue.GetLocationIndexUsingCell(&(spring_iterator.rGetCellB())),
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

    void TestGetLocationOfCellAndGetNodeCorrespondingToCell() throw (Exception)
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node. Get each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create the tissue
        MeshBasedTissue<2> tissue(mesh, cells);

        // Loop over nodes
        for (AbstractTissue<2>::Iterator cell_iter=tissue.Begin();
             cell_iter!=tissue.End();
             ++cell_iter)
        {
            // Record node location
            c_vector<double,2> node_location = tissue.GetLocationOfCell(&(*cell_iter));

            // Test GetLocationOfCell()
            TS_ASSERT_DELTA(node_location[0], tissue.GetLocationOfCell(&(*cell_iter))[0], 1e-9);
            TS_ASSERT_DELTA(node_location[1], tissue.GetLocationOfCell(&(*cell_iter))[1], 1e-9);
        }
    }

    // At the moment the tissue cannot be properly archived since the mesh cannot be. This test
    // just checks that the cells are correctly archived.
    void TestArchivingMeshBasedTissue() throw (Exception)
    {
        OutputFileHandler handler("archive",false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "tissue.arch";

        // Archive a tissue
        {
            // Need to set up time
            unsigned num_steps=10;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);

            // Create a simple mesh
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
            MutableMesh<2,2> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            // Set up cells, one for each node. Get each a birth time of -node_index,
            // so the age = node_index
            std::vector<TissueCell> cells;
            FixedCellCycleModelCellsGenerator<2> cells_generator;
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
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);
            p_simulation_time->IncrementTimeOneStep();

            MeshBasedTissue<2>* p_tissue;

            // Restore the tissue
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> *p_simulation_time;

            // The following line is required because the loading of a tissue
            // is usually called by the method TissueSimulation::Load()
            MeshArchiveInfo::meshPathname = "mesh/test/data/square_4_elements";

            input_arch >> p_tissue;

            // Cells have been given birth times of 0, -1, -2, -3, -4.
            // this checks that individual cells and their models are archived.
            unsigned counter = 0u;
            for (AbstractTissue<2>::Iterator cell_iter=p_tissue->Begin();
                 cell_iter!=p_tissue->End();
                 ++cell_iter)
            {
                TS_ASSERT_DELTA(cell_iter->GetAge(),(double)(counter),1e-7);
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

            // This won't pass because of the mesh not being archived
            // TS_ASSERT_EQUALS(tissue.rGetMesh().GetNumNodes(),5u);

            delete p_tissue;
        }
    }

    void TestSpringMarking()
    {
        // Create a small tissue
        std::vector<Node<2> *> nodes;
        nodes.push_back(new Node<2>(0, false, 0, 0.5));
        nodes.push_back(new Node<2>(1, false, 1, 0));
        nodes.push_back(new Node<2>(2, false, 1, 1));
        nodes.push_back(new Node<2>(3, false, 2, 0.5));
        nodes.push_back(new Node<2>(4, false, 2, 1.5));

        MutableMesh<2,2> mesh(nodes);

        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
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
        AbstractTissue<2>::Iterator it=tissue.Begin();
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
        std::vector<Node<2> *> nodes;
        nodes.push_back(new Node<2>(0, false, 0, 0.5));
        nodes.push_back(new Node<2>(1, false, 1, 0));
        nodes.push_back(new Node<2>(2, false, 1, 1));
        nodes.push_back(new Node<2>(3, false, 2, 0.5));
        nodes.push_back(new Node<2>(4, false, 2, 1.5));

        MutableMesh<2,2> mesh(nodes);

        std::vector<TissueCell> cells;
        FixedCellCycleModelCellsGenerator<2> cells_generator;
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
