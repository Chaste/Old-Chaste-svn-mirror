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
#ifndef TESTNODEBASEDTISSUE_HPP_
#define TESTNODEBASEDTISSUE_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "NodeBasedTissue.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "TrianglesMeshReader.hpp"
#include "TetrahedralMesh.hpp"
#include "AbstractCellBasedTestSuite.hpp"


class TestNodeBasedTissue : public AbstractCellBasedTestSuite
{
private:

    template<unsigned DIM>
    std::vector<TissueCell> SetUpCells(TetrahedralMesh<DIM,DIM>* pMesh)
    {
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<pMesh->GetNumNodes(); i++)
        {
            TissueCell cell(STEM, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
            double birth_time = 0.0 - i;
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        return cells;
    }

    template<unsigned DIM>
    void TestSimpleNodeBasedTissue(std::string meshFilename)
    {
        // Create a simple mesh
        TrianglesMeshReader<DIM,DIM> mesh_reader(meshFilename);
        TetrahedralMesh<DIM,DIM> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node. Give each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells = SetUpCells(&mesh);

        // Create the tissue
        NodeBasedTissue<DIM> node_based_tissue(mesh, cells);

        TS_ASSERT_EQUALS(node_based_tissue.rGetNodes().size(), mesh.GetNumNodes());
        TS_ASSERT_EQUALS(node_based_tissue.rGetCells().size(), cells.size());

        unsigned counter = 0;
        for (typename AbstractTissue<DIM>::Iterator cell_iter = node_based_tissue.Begin();
             cell_iter != node_based_tissue.End();
             ++cell_iter)
        {
            // Test operator* and that cells are in sync
            TS_ASSERT_EQUALS(node_based_tissue.GetLocationIndexUsingCell(*cell_iter), counter);

            // Test operator-> and that cells are in sync
            TS_ASSERT_DELTA(cell_iter->GetAge(), (double)counter, 1e-12);

            counter++;
        }

        TS_ASSERT_EQUALS(counter, node_based_tissue.GetNumRealCells());
    }


public:

    // Test construction, accessors and Iterator
    void TestNodeBasedTissue1d2d3d() throw(Exception)
    {
        TestSimpleNodeBasedTissue<1>("mesh/test/data/1D_0_to_1_10_elements");
        TestSimpleNodeBasedTissue<2>("mesh/test/data/square_4_elements");
        TestSimpleNodeBasedTissue<3>("mesh/test/data/cube_136_elements");
    }

    void TestOtherNodeBasedTissueConstructor()
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node. Give each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells = SetUpCells(&mesh);

        // Get a std::vector of nodes from the mesh
        std::vector<Node<2>* > nodes;

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            Node<2>* p_node = new Node<2>(*(mesh.GetNode(i)));
            nodes.push_back(p_node);
        }

        // Create the tissue
        NodeBasedTissue<2> node_based_tissue(nodes, cells);

        TS_ASSERT_EQUALS(node_based_tissue.rGetNodes().size(), mesh.GetNumNodes());
        TS_ASSERT_EQUALS(node_based_tissue.rGetNodes().size(), nodes.size());
        TS_ASSERT_EQUALS(node_based_tissue.rGetCells().size(), cells.size());

        // For coverage, test that tissue constructor with 3rd argument locationIndices throws
        // an exception when the size of locationIndices does not equal the number of cells
        std::vector<unsigned> location_indices;
        location_indices.push_back(0);
        location_indices.push_back(1);
        location_indices.push_back(2);

        TS_ASSERT_THROWS_THIS(NodeBasedTissue<2> node_based_tissue(nodes, cells, location_indices),"There is not a one-one correspondence between cells and location indices");
    }

    void TestValidateNodeBasedTissue()
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node apart from one.
        // Give each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<mesh.GetNumNodes()-1; i++)
        {
            AbstractCellCycleModel* p_cell_cycle_model = new FixedDurationGenerationBasedCellCycleModel();
            TissueCell cell(STEM, HEALTHY, p_cell_cycle_model);
            double birth_time = 0.0 - i;
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Get a std::vector of nodes from the mesh
        std::vector<Node<2>* > nodes;

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            Node<2>* p_node = new Node<2>(*(mesh.GetNode(i)));
            nodes.push_back(p_node);
        }
        // Fails as no cell corresponding to node 4
        TS_ASSERT_THROWS_THIS(NodeBasedTissue<2> tissue(nodes, cells),
                "Node 4 does not appear to have a cell associated with it");

        // Add another cell
        AbstractCellCycleModel* p_cell_cycle_model = new FixedDurationGenerationBasedCellCycleModel();
        TissueCell cell(STEM, HEALTHY, p_cell_cycle_model);
        double birth_time = -4.0;
        cell.SetBirthTime(birth_time);
        cells.push_back(cell);

        NodeBasedTissue<2> tissue(nodes, cells);

        // throws exception as update hasn't been called so no node pairs set up yet
        TS_ASSERT_THROWS_THIS(tissue.rGetNodePairs(),
                "No node pairs set up, rGetNodePairs probably called before Update");

        // throws exception as the cut-off length hasn't been set and has its default value of DBL_MAX
        TS_ASSERT_THROWS_THIS(tissue.Update(),
                "NodeBasedTissue cannot create boxes if the cut-off length has not been set - Call UseCutoffPoint() on the force law, or SetMechanicsCutOffLength on TissueConfig");

        TissueConfig::Instance()->SetMechanicsCutOffLength(1.2);
        tissue.Update();

        std::set< std::pair<Node<2>*, Node<2>* > >& r_node_pairs = tissue.rGetNodePairs();
        r_node_pairs.clear();
    }

    void TestAddCell()
    {
        // Create two nodes
        ChastePoint<2> point0;
        point0.rGetLocation()[0] = 0.0;
        point0.rGetLocation()[1] = 0.0;
        Node<2>* p_node0 = new Node<2>(0, point0, false);

        ChastePoint<2> point1;
        point1.rGetLocation()[0] = 1.0;
        point1.rGetLocation()[1] = 1.0;
        Node<2>* p_node1 = new Node<2>(1, point1, false);

        std::vector<Node<2>* > nodes;
        nodes.push_back(p_node0);
        nodes.push_back(p_node1);

        // Create two cells
        TissueCell cell0(STEM, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
        cell0.SetBirthTime(-1);

        TissueCell cell1(STEM, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
        cell1.SetBirthTime(-1);

        std::vector<TissueCell> cells;
        cells.push_back(cell0);
        cells.push_back(cell1);

        // Create a tissue
        NodeBasedTissue<2> node_based_tissue(nodes, cells);

        // Create a new cell, DON'T set the node index, set birth time=-1
        TissueCell cell2(STEM, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
        cell2.SetBirthTime(-1);

        c_vector<double,2> cell2_location;
        cell2_location[0] = 2.0;
        cell2_location[1] = 2.0;

        node_based_tissue.AddCell(cell2, cell2_location);
    }

    void TestSetNodeAndAddCell()
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node. Give each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells = SetUpCells(&mesh);

        // Create a tissue
        NodeBasedTissue<2> node_based_tissue(mesh, cells);

        // Test SetNode() by moving node 0 by a small amount

        AbstractTissue<2>::Iterator cell_iter = node_based_tissue.Begin();
        c_vector<double,2> new_location = node_based_tissue.GetLocationOfCellCentre(*cell_iter);
        new_location[0] += 1e-2;
        new_location[1] += 1e-2;
        ChastePoint<2> new_location_point(new_location);
        node_based_tissue.SetNode(node_based_tissue.GetLocationIndexUsingCell(*cell_iter), new_location_point);

        TS_ASSERT_DELTA(node_based_tissue.GetNode(0)->rGetLocation()[0], new_location[0], 1e-12);
        TS_ASSERT_DELTA(node_based_tissue.GetNode(0)->rGetLocation()[1], new_location[1], 1e-12);

        // Test AddNode

        ChastePoint<2> new_point;
        new_point.rGetLocation()[0] = 1.71;
        new_point.rGetLocation()[1] = 1.72;

        unsigned num_nodes = node_based_tissue.GetNumNodes();

        Node<2>* p_node = new Node<2>(num_nodes, new_point, false);
        unsigned new_node_index = node_based_tissue.AddNode(p_node);

        TS_ASSERT_EQUALS(new_node_index, num_nodes);
        TS_ASSERT_DELTA(node_based_tissue.GetNode(num_nodes)->rGetLocation()[0], 1.71, 1e-12);
        TS_ASSERT_DELTA(node_based_tissue.GetNode(num_nodes)->rGetLocation()[1], 1.72, 1e-12);

        // Tidy up
        node_based_tissue.mNodes.pop_back();
        delete p_node;

        // Remove a cell so as to populate mDeletedNodeIndices (for coverage)
        node_based_tissue.rGetCellUsingLocationIndex(0).Kill();
        node_based_tissue.RemoveDeadCells();

        // Test AddNode again

        ChastePoint<2> new_point2;
        new_point2.rGetLocation()[0] = 0.51;
        new_point2.rGetLocation()[1] = 0.52;

        num_nodes = node_based_tissue.GetNumNodes();
        Node<2>* p_node2 = new Node<2>(num_nodes, new_point2, false);
        new_node_index = node_based_tissue.AddNode(p_node2);

        TS_ASSERT_EQUALS(new_node_index, 0u);
        TS_ASSERT_DELTA(node_based_tissue.GetNode(0)->rGetLocation()[0], 0.51, 1e-12);
        TS_ASSERT_DELTA(node_based_tissue.GetNode(0)->rGetLocation()[1], 0.52, 1e-12);

        // Test AddCell

        unsigned old_num_nodes = node_based_tissue.rGetNodes().size();
        unsigned old_num_cells = node_based_tissue.rGetCells().size();

        // Create a new cell, DON'T set the node index, set birth time=-1
        TissueCell cell(STEM, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
        cell.SetBirthTime(-1);

        c_vector<double,2> new_cell_location;
        new_cell_location[0] = 2.0;
        new_cell_location[1] = 2.0;

        node_based_tissue.AddCell(cell, new_cell_location);

        // Tissue should have updated nodes and cells
        TS_ASSERT_EQUALS(node_based_tissue.GetNumNodes(), old_num_nodes+1);
        TS_ASSERT_EQUALS(node_based_tissue.rGetCells().size(), old_num_cells+1);
        TS_ASSERT_EQUALS(node_based_tissue.GetNumRealCells(), old_num_nodes);

        // Check the location of the new node
        TS_ASSERT_DELTA(node_based_tissue.GetNode(old_num_nodes)->rGetLocation()[0], 2.0, 1e-12);
        TS_ASSERT_DELTA(node_based_tissue.GetNode(old_num_nodes)->rGetLocation()[1], 2.0, 1e-12);

        // Check the index of the new cell
        TissueCell& new_cell = node_based_tissue.rGetCells().back();
        TS_ASSERT_EQUALS(node_based_tissue.GetLocationIndexUsingCell(new_cell), old_num_nodes);
    }

    void TestRemoveDeadCellsAndUpdate()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node. Give each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells = SetUpCells(&mesh);

        cells[27].StartApoptosis();

        // Create a tissue
        NodeBasedTissue<2> node_based_tissue(mesh, cells);

        // Test we have the right numbers of nodes and cells
        TS_ASSERT_EQUALS(node_based_tissue.GetNumNodes(), 81u);
        TS_ASSERT_EQUALS(node_based_tissue.GetNumRealCells(), 81u);

        p_simulation_time->IncrementTimeOneStep();

        unsigned num_removed = node_based_tissue.RemoveDeadCells();

        TissueConfig::Instance()->SetMechanicsCutOffLength(1.2);
        node_based_tissue.Update(true);

        // Test that one cell has been removed
        TS_ASSERT_EQUALS(num_removed, 1u);
        TS_ASSERT_EQUALS(node_based_tissue.GetNumRealCells(), 80u);

        // Test that one node has been removed
        TS_ASSERT_EQUALS(node_based_tissue.GetNumNodes(), 80u);

        // Test that each cell's node index has been correctly updated
        unsigned index = 0;
        for (AbstractTissue<2>::Iterator cell_iter = node_based_tissue.Begin();
             cell_iter != node_based_tissue.End();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS(node_based_tissue.GetLocationIndexUsingCell(*cell_iter), index);
            index++;
        }
    }

    void TestAddAndRemoveAndAddWithOutRemovingDeletedNodes()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node. Give each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells = SetUpCells(&mesh);

        // Make one cell start apoptosis
        cells[27].StartApoptosis();

        // Create a tissue
        NodeBasedTissue<2> node_based_tissue(mesh, cells);

        // Test we have the right numbers of nodes and cells
        TS_ASSERT_EQUALS(node_based_tissue.GetNumNodes(), 81u);
        TS_ASSERT_EQUALS(node_based_tissue.GetNumRealCells(), 81u);

        // Add a cell to the tissue
        TissueCell new_cell(STEM, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
        new_cell.SetBirthTime(0);
        c_vector<double,2> new_location;
        new_location[0] = 0.3433453454443;
        new_location[1] = 0.3435346344234;
        node_based_tissue.AddCell(new_cell, new_location);

        // Test that the numbers of nodes and cells has been updated
        TS_ASSERT_EQUALS(node_based_tissue.GetNumNodes(), 82u);
        TS_ASSERT_EQUALS(node_based_tissue.GetNumRealCells(), 82u);

        p_simulation_time->IncrementTimeOneStep();

        // Test that the apoptotic cell has been removed
        unsigned num_removed = node_based_tissue.RemoveDeadCells();
        TissueConfig::Instance()->SetMechanicsCutOffLength(1.2);
        node_based_tissue.Update();

        TS_ASSERT_EQUALS(num_removed, 1u);
        TS_ASSERT_EQUALS(node_based_tissue.GetNumNodes(), 81u);
        TS_ASSERT_EQUALS(node_based_tissue.GetNumRealCells(), 81u);

        // Add another cell to the tissue
        TissueCell new_cell2(STEM, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
        new_cell2.SetBirthTime(0);

        c_vector<double,2> new_location2;
        new_location2[0] = 0.6433453454443;
        new_location2[1] = 0.6435346344234;
        node_based_tissue.AddCell(new_cell2, new_location2);

        // Test that the numbers of nodes and cells has been updated
        TS_ASSERT_EQUALS(node_based_tissue.GetNumNodes(), 82u);
        TS_ASSERT_EQUALS(node_based_tissue.GetNumRealCells(), 82u);
    }

    void TestSettingCellAncestors() throw (Exception)
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node. Give each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells = SetUpCells(&mesh);

        // Create a tissue
        NodeBasedTissue<2> node_based_tissue(mesh, cells);

        // Test that the tissue makes all cells fix the node index as ancestor
        node_based_tissue.SetCellAncestorsToLocationIndices();

        unsigned counter = 0;

        for (AbstractTissue<2>::Iterator cell_iter = node_based_tissue.Begin();
             cell_iter != node_based_tissue.End();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS(cell_iter->GetAncestor(), node_based_tissue.GetLocationIndexUsingCell(*cell_iter));
            counter ++;
        }
        TS_ASSERT_EQUALS(counter, 5u);

        // Test that we can recover remaining number of ancestors
        std::set<unsigned> remaining_ancestors = node_based_tissue.GetCellAncestors();
        TS_ASSERT_EQUALS(remaining_ancestors.size(), 5u);

        // Test that the set correctly represents a monoclonal population
        for (AbstractTissue<2>::Iterator cell_iter=node_based_tissue.Begin();
             cell_iter!=node_based_tissue.End();
             ++cell_iter)
        {
            // Set all cells to have the same ancestor...
            cell_iter->SetAncestor(1u);
        }
        remaining_ancestors = node_based_tissue.GetCellAncestors();
        TS_ASSERT_EQUALS(remaining_ancestors.size(), 1u);
    }

    void TestGetLocationOfCellCentre() throw (Exception)
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node. Give each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells = SetUpCells(&mesh);

        // Create a tissue
        NodeBasedTissue<2> node_based_tissue(mesh, cells);

        // Loop over nodes
        for (AbstractTissue<2>::Iterator cell_iter = node_based_tissue.Begin();
             cell_iter != node_based_tissue.End();
             ++cell_iter)
        {
            // Record node location
            c_vector<double, 2> node_location = node_based_tissue.GetLocationOfCellCentre(*cell_iter);

            // Test GetLocationOfCellCentre()
            TS_ASSERT_DELTA(node_location[0], node_based_tissue.GetLocationOfCellCentre(*cell_iter)[0], 1e-9);
            TS_ASSERT_DELTA(node_location[1], node_based_tissue.GetLocationOfCellCentre(*cell_iter)[1], 1e-9);
        }
    }

    void TestNodeBasedTissueOutputWriters()
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node. Give each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells = SetUpCells(&mesh);

        // Create a tissue
        NodeBasedTissue<2> node_based_tissue(mesh, cells);

        // For coverage of WriteResultsToFiles()
        node_based_tissue.rGetCellUsingLocationIndex(0).SetCellProliferativeType(TRANSIT);
        node_based_tissue.rGetCellUsingLocationIndex(0).SetMutationState(LABELLED);
        node_based_tissue.rGetCellUsingLocationIndex(1).SetCellProliferativeType(DIFFERENTIATED);
        node_based_tissue.rGetCellUsingLocationIndex(1).SetMutationState(APC_ONE_HIT);
        node_based_tissue.rGetCellUsingLocationIndex(2).SetMutationState(APC_TWO_HIT);
        node_based_tissue.rGetCellUsingLocationIndex(3).SetMutationState(BETA_CATENIN_ONE_HIT);
        node_based_tissue.rGetCellUsingLocationIndex(4).SetCellProliferativeType(APOPTOTIC);
        node_based_tissue.rGetCellUsingLocationIndex(4).StartApoptosis();
        node_based_tissue.SetCellAncestorsToLocationIndices();

        std::string output_directory = "TestNodeBasedTissueWriters";
        OutputFileHandler output_file_handler(output_directory, false);

        TissueConfig::Instance()->SetOutputCellMutationStates(true);
        TissueConfig::Instance()->SetOutputCellProliferativeTypes(true);
        TissueConfig::Instance()->SetOutputCellCyclePhases(true);
        TissueConfig::Instance()->SetOutputCellAncestors(true);
        TissueConfig::Instance()->SetOutputCellAges(true);

        TS_ASSERT_THROWS_NOTHING(node_based_tissue.CreateOutputFiles(output_directory, false));

        node_based_tissue.WriteResultsToFiles();

        TS_ASSERT_THROWS_NOTHING(node_based_tissue.CloseOutputFiles());

        // Compare output with saved files of what they should look like
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.viznodes     cancer/test/data/TestNodeBasedTissueWriters/results.viznodes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizcelltypes     cancer/test/data/TestNodeBasedTissueWriters/results.vizcelltypes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizancestors     cancer/test/data/TestNodeBasedTissueWriters/results.vizancestors").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellmutationstates.dat     cancer/test/data/TestNodeBasedTissueWriters/cellmutationstates.dat").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellages.dat     cancer/test/data/TestNodeBasedTissueWriters/cellages.dat").c_str()), 0);

        // Test the GetCellMutationStateCount function
        std::vector<unsigned> cell_mutation_states = node_based_tissue.rGetCellMutationStateCount();
        TS_ASSERT_EQUALS(cell_mutation_states.size(), 5u);
         TS_ASSERT_EQUALS(cell_mutation_states[0], 1u);
        for (unsigned i=1; i<cell_mutation_states.size(); i++)
        {
            TS_ASSERT_EQUALS(cell_mutation_states[i], 1u);
        }

        // Test the GetCellProliferativeTypeCount function
        std::vector<unsigned> cell_types = node_based_tissue.rGetCellProliferativeTypeCount();
        TS_ASSERT_EQUALS(cell_types.size(), 4u);
        TS_ASSERT_EQUALS(cell_types[0], 2u);
        TS_ASSERT_EQUALS(cell_types[1], 1u);
        TS_ASSERT_EQUALS(cell_types[2], 1u);
        TS_ASSERT_EQUALS(cell_types[3], 1u);

        // For coverage
        TS_ASSERT_THROWS_NOTHING(node_based_tissue.WriteResultsToFiles());
    }

    void TestWritingCellCyclePhases()
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<TissueCell> cells;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            TissueCell cell(STEM, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
            double birth_time;
            if (i==1)
            {
                birth_time = -0.5;
            }
            else if (i==2)
            {
                birth_time = -1.5;
            }
            else if (i==3)
            {
                birth_time = -15.5;
            }
            else
            {
                birth_time = -23.5;
            }
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        cells[0].SetCellProliferativeType(DIFFERENTIATED);

        // Create a tissue
        NodeBasedTissue<2> node_based_tissue(mesh, cells);

        // Cells have been given birth times of 0, -1, -2, -3, -4.
        // loop over them to run to time 0.0;
        for (AbstractTissue<2>::Iterator cell_iter = node_based_tissue.Begin();
             cell_iter != node_based_tissue.End();
             ++cell_iter)
        {
            cell_iter->ReadyToDivide();
        }

        std::string output_directory = "TestWritingCellCyclePhases";
        OutputFileHandler output_file_handler(output_directory, false);

        TissueConfig::Instance()->SetOutputCellCyclePhases(true);

        node_based_tissue.CreateOutputFiles(output_directory, false);
        node_based_tissue.WriteResultsToFiles();
        node_based_tissue.CloseOutputFiles();

        // Test the rGetCellCyclePhaseCount function
        std::vector<unsigned> cell_cycle_phases = node_based_tissue.rGetCellCyclePhaseCount();
        TS_ASSERT_EQUALS(cell_cycle_phases[0], 1u);
        TS_ASSERT_EQUALS(cell_cycle_phases[1], 3u);
        TS_ASSERT_EQUALS(cell_cycle_phases[2], 0u);
        TS_ASSERT_EQUALS(cell_cycle_phases[3], 0u);
        TS_ASSERT_EQUALS(cell_cycle_phases[4], 1u);
    }

    void TestArchivingTissue() throw (Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "NodeBasedTissue.arch";

        // Archive a simple tissue
        {
            // Need to set up time
            unsigned num_steps = 10;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);

            // Create a simple mesh
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
            TetrahedralMesh<2,2> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            // Set up cells, one for each node. Give each a birth time of -node_index,
            // so the age = node_index
            std::vector<TissueCell> cells = SetUpCells(&mesh);

            // Create a tissue
            NodeBasedTissue<2>* const p_tissue = new NodeBasedTissue<2>(mesh, cells);

            // Cells have been given birth times of 0, -1, -2, -3, -4.
            // loop over them to run to time 0.0;
            for (AbstractTissue<2>::Iterator cell_iter = p_tissue->Begin();
                cell_iter != p_tissue->End();
                ++cell_iter)
            {
                cell_iter->ReadyToDivide();
            }

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Write the tissue to the archive
            output_arch << static_cast<const SimulationTime&> (*p_simulation_time);
            output_arch << p_tissue;
            SimulationTime::Destroy();
            delete p_tissue;
        }

        // Restore simple tissue
        {
            // Need to set up time
            unsigned num_steps = 10;

            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);
            p_simulation_time->IncrementTimeOneStep();

            NodeBasedTissue<2>* p_tissue;

            // Restore the tissue
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> *p_simulation_time;
            input_arch >> p_tissue;

            // Cells have been given birth times of 0, -1, -2, -3, -4.
            // this checks that individual cells and their models are archived.
            unsigned counter = 0;

            for (AbstractTissue<2>::Iterator cell_iter = p_tissue->Begin();
                 cell_iter != p_tissue->End();
                 ++cell_iter)
            {
                TS_ASSERT_DELTA(cell_iter->GetAge(), (double)(counter), 1e-7);
                counter++;
            }

            // Check the simulation time has been restored (through the cell)
            TS_ASSERT_EQUALS(p_simulation_time->GetTime(), 0.0);

            // Check the tissue has been restored
            TS_ASSERT_EQUALS(p_tissue->rGetCells().size(), 5u);

            // Check number of nodes
            std::vector<Node<2>* > nodes = p_tissue->rGetNodes();
            TS_ASSERT_EQUALS(nodes.size(), 5u);

            // Check some node positions
            TS_ASSERT_EQUALS(nodes[3]->GetIndex(), 3u);
            TS_ASSERT_EQUALS(nodes[4]->GetIndex(), 4u);

            TS_ASSERT_DELTA(nodes[3]->rGetLocation()[0], 0.0, 1e-9);
            TS_ASSERT_DELTA(nodes[3]->rGetLocation()[1], 1.0, 1e-9);
            TS_ASSERT_DELTA(nodes[4]->rGetLocation()[0], 0.5, 1e-9);
            TS_ASSERT_DELTA(nodes[4]->rGetLocation()[1], 0.5, 1e-9);

            TS_ASSERT_EQUALS(nodes[3]->IsBoundaryNode(), true);
            TS_ASSERT_EQUALS(nodes[4]->IsBoundaryNode(), false);

            delete p_tissue;
        }
    }



};

#endif /*TESTNODEBASEDTISSUE_HPP_*/
