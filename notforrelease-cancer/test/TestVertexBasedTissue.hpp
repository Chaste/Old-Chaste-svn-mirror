/*

Copyright (C) University of Oxford, 2008

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

#ifndef TESTVERTEXBASEDTISSUE_HPP_
#define TESTVERTEXBASEDTISSUE_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "VertexBasedTissue.hpp"
#include "FixedCellCycleModel.hpp"
#include "AbstractCancerTestSuite.hpp"


class TestVertexBasedTissue : public AbstractCancerTestSuite
{
private:

    /**
     * Set up cells, one for each VertexElement.
     * Give each cell a birth time of -elem_index,
     * so its age is elem_index.
     */
    template<unsigned DIM>
    std::vector<TissueCell> SetUpCells(VertexMesh<DIM,DIM>& rMesh)
    {
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<rMesh.GetNumElements(); i++)
        {
            TissueCell cell(DIFFERENTIATED, HEALTHY, new FixedCellCycleModel());
            double birth_time = 0.0 - i;
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        return cells;
    }

public:

    // Test construction, accessors and iterator
    void TestCreateSmallVertexBasedTissue() throw(Exception)
    {
        // Create a simple 2D VertexMesh
        VertexMesh<2,2> mesh(5,3); // columns then rows

        // Set up cells
        std::vector<TissueCell> cells = SetUpCells(mesh);

        // Create tissue
        VertexBasedTissue<2> tissue(mesh, cells);

        // Test we have the correct number of cells and elements
        TS_ASSERT_EQUALS(tissue.GetNumElements(), mesh.GetNumElements());
        TS_ASSERT_EQUALS(tissue.rGetCells().size(), cells.size());

        unsigned counter = 0;

        // Test VertexBasedTissue::Iterator
        for (VertexBasedTissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            // Test operator* and that cells are in sync
            TS_ASSERT_EQUALS(tissue.GetElementCorrespondingToCell(&(*cell_iter))->GetIndex(), counter);

            // Test operator-> and that cells are in sync
            TS_ASSERT_DELTA(cell_iter->GetAge(), (double)counter, 1e-12);

            TS_ASSERT_EQUALS(counter, mesh.GetElement(counter)->GetIndex());

            counter++;
        }

        // Test we have gone through all cells in the for loop
        TS_ASSERT_EQUALS(counter, tissue.GetNumRealCells());

        // Test GetNumNodes() method
        TS_ASSERT_EQUALS(tissue.GetNumNodes(), mesh.GetNumNodes());
    }

    void TestValidateVertexBasedTissue()
    {
        // Create a simple vertex-based mesh
        VertexMesh<2,2> mesh(4,6); // columns then rows

        // Set up cells, one for each element.
        // Get each a birth time of -element_index, so the age = element_index.
        std::vector<TissueCell> cells;
        std::vector<unsigned> cell_location_indices;
        for (unsigned i=0; i<mesh.GetNumElements()-1; i++)
        {
            AbstractCellCycleModel* p_cell_cycle_model = new FixedCellCycleModel();
            TissueCell cell(STEM, HEALTHY, p_cell_cycle_model);
            double birth_time = 0.0 - i;
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
            
            cell_location_indices.push_back(i);
        }

        TS_ASSERT_THROWS_ANYTHING(VertexBasedTissue<2> tissue(mesh, cells));

        AbstractCellCycleModel* p_cell_cycle_model = new FixedCellCycleModel();
        TissueCell cell(STEM, HEALTHY, p_cell_cycle_model);
        double birth_time = 0.0 - mesh.GetNumElements()-1;
        cell.SetBirthTime(birth_time);
        cells.push_back(cell);        
        cell_location_indices.push_back(mesh.GetNumElements()-1);

        TS_ASSERT_THROWS_NOTHING(VertexBasedTissue<2> tissue(mesh, cells));
    }

    void TestUpdateWithoutBirthOrDeath()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);

        // Create a simple vertex-based mesh
        VertexMesh<2,2> mesh(4,6); // columns then rows

        // Set up cells
        std::vector<TissueCell> cells = SetUpCells(mesh);

        // Create tissue
        VertexBasedTissue<2> tissue(mesh, cells);

        /// \todo Coverage (can be removed once test below is completed - see #853)
        unsigned num_cells_removed = tissue.RemoveDeadCells();
        TS_ASSERT_EQUALS(num_cells_removed, 0u);

        p_simulation_time->IncrementTimeOneStep();

        TS_ASSERT_THROWS_NOTHING(tissue.Update());
    }

 
    void TestAddCell()
    {
        // Make four nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 2.0, -1.0));
        nodes.push_back(new Node<2>(1, false, 2.0, 1.0));
        nodes.push_back(new Node<2>(2, false, -2.0, 1.0));
        nodes.push_back(new Node<2>(3, false, -2.0, -1.0));
        nodes.push_back(new Node<2>(4, false, 0.0, 2.0));

        // Make a rectangular element out of nodes 0,1,2,3
        std::vector<Node<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[0]);
        nodes_elem_1.push_back(nodes[1]);
        nodes_elem_1.push_back(nodes[2]);
        nodes_elem_1.push_back(nodes[3]);

        // Make a triangular element out of nodes 1,4,2
        std::vector<Node<2>*> nodes_elem_2;
        nodes_elem_2.push_back(nodes[1]);
        nodes_elem_2.push_back(nodes[4]);
        nodes_elem_2.push_back(nodes[2]);

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_1));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_2));

        // Make a vertex mesh
        VertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 5u);

        // Create cells
        std::vector<TissueCell> cells = SetUpCells(vertex_mesh);

        // Create tissue
        VertexBasedTissue<2> tissue(vertex_mesh, cells);

        unsigned old_num_nodes = vertex_mesh.GetNumNodes();
        unsigned old_num_elements = vertex_mesh.GetNumElements();
        unsigned old_num_cells = tissue.rGetCells().size();

        // Add new cell by dividing element 0 along short axis
       
        c_vector<double,2> new_cell_location;
        new_cell_location[0] = 2.0;
        new_cell_location[1] = 2.0;

        TissueCell cell0 = tissue.rGetCellUsingLocationIndex(0);

        /// \todo Second input argument is not needed - tidy this up (#852)        
        TissueCell new_cell(STEM, HEALTHY, new FixedCellCycleModel());
        new_cell.SetBirthTime(-1);
        
        TissueCell* p_new_cell = tissue.AddCell(new_cell, new_cell_location, &cell0);

        // Check that the new cell was successfully added to the tissue
        TS_ASSERT_EQUALS(tissue.GetNumNodes(), old_num_nodes+2);
        TS_ASSERT_EQUALS(tissue.GetNumElements(), old_num_elements+1);
        TS_ASSERT_EQUALS(tissue.rGetCells().size(), old_num_cells+1);
        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), old_num_elements+1);

        // Check the location of the new nodes
        TS_ASSERT_DELTA(tissue.GetNode(old_num_nodes)->rGetLocation()[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(tissue.GetNode(old_num_nodes)->rGetLocation()[1], 1.0, 1e-12);

        TS_ASSERT_DELTA(tissue.GetNode(old_num_nodes+1)->rGetLocation()[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(tissue.GetNode(old_num_nodes+1)->rGetLocation()[1], -1.0, 1e-12);

        // Now test the nodes in each element
        TS_ASSERT_EQUALS(tissue.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(tissue.GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(tissue.GetElement(2)->GetNumNodes(), 4u);

        TS_ASSERT_EQUALS(tissue.GetElement(0)->GetNodeGlobalIndex(0), 5u);
        TS_ASSERT_EQUALS(tissue.GetElement(0)->GetNodeGlobalIndex(1), 2u);
        TS_ASSERT_EQUALS(tissue.GetElement(0)->GetNodeGlobalIndex(2), 3u);
        TS_ASSERT_EQUALS(tissue.GetElement(0)->GetNodeGlobalIndex(3), 6u);

        TS_ASSERT_EQUALS(tissue.GetElement(1)->GetNodeGlobalIndex(0), 1u);
        TS_ASSERT_EQUALS(tissue.GetElement(1)->GetNodeGlobalIndex(1), 4u);
        TS_ASSERT_EQUALS(tissue.GetElement(1)->GetNodeGlobalIndex(2), 2u);
        TS_ASSERT_EQUALS(tissue.GetElement(1)->GetNodeGlobalIndex(3), 5u);

        TS_ASSERT_EQUALS(tissue.GetElement(2)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(tissue.GetElement(2)->GetNodeGlobalIndex(1), 1u);
        TS_ASSERT_EQUALS(tissue.GetElement(2)->GetNodeGlobalIndex(2), 5u);
        TS_ASSERT_EQUALS(tissue.GetElement(2)->GetNodeGlobalIndex(3), 6u);
        
        // Test ownership of the new nodes
        std::set<unsigned> expected_elements_containing_node_5;
        expected_elements_containing_node_5.insert(0);
        expected_elements_containing_node_5.insert(1);
        expected_elements_containing_node_5.insert(2);
        
        TS_ASSERT_EQUALS(tissue.GetNode(5)->rGetContainingElementIndices(), expected_elements_containing_node_5);
        
        std::set<unsigned> expected_elements_containing_node_6;
        expected_elements_containing_node_6.insert(0);
        expected_elements_containing_node_6.insert(2);
        
        TS_ASSERT_EQUALS(tissue.GetNode(6)->rGetContainingElementIndices(), expected_elements_containing_node_6);

        // Check the index of the new cell
        TS_ASSERT_EQUALS(tissue.GetElementCorrespondingToCell(p_new_cell)->GetIndex(), old_num_elements);
    }


    /// \todo This test currently fails, since the method RemoveDeadCells() does not yet
    // delete the elements/nodes associated with dead cells (see #853)
    void DONTTestRemoveDeadCellsAndUpdate()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);

        // Create a simple vertex-based mesh
        VertexMesh<2,2> mesh(4,6); // columns then rows

        // Set up cells
        std::vector<TissueCell> cells = SetUpCells(mesh);
        cells[5].StartApoptosis();

        // Create tissue
        VertexBasedTissue<2> tissue(mesh, cells);

        TS_ASSERT_EQUALS(mesh.GetNumElements(), 24u);
        TS_ASSERT_EQUALS(tissue.rGetCells().size(), 24u);
        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), 24u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 68u);

        p_simulation_time->IncrementTimeOneStep();

        // Remove dead cells
        unsigned num_cells_removed = tissue.RemoveDeadCells();

        /// \todo Currently RemoveDeadCells() does nothing, and is only
        //        in a test for coverage. Cell death will be implemented
        //        in #853.
        TS_ASSERT_EQUALS(num_cells_removed, 1u);

        // We should now have one less real cell, since one cell has been
        // marked as dead, so is skipped by the tissue iterator
        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), 23u);

        /// \todo Need some more tests here, on the new number of elements/nodes

        TS_ASSERT_EQUALS(tissue.rGetCells().size(), cells.size()); // the tissue now copies cells

        tissue.Update();

        // Finally, check the cells' element indices have updated

        // We expect the cell element indices to be {0,11,...,23}
        std::set<unsigned> expected_elem_indices;
        for (unsigned i=0; i<tissue.GetNumRealCells(); i++)
        {
            expected_elem_indices.insert(i);
        }

        // Get actual cell element indices
        std::set<unsigned> element_indices;

        for (AbstractTissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {
            // Record element index corresponding to cell
            unsigned element_index = tissue.GetElementCorrespondingToCell(&(*cell_iter))->GetIndex();
            element_indices.insert(element_index);
        }

        TS_ASSERT_EQUALS(element_indices, expected_elem_indices);
    }

    void TestVertexBasedTissueOutputWriters()
    {
        // Create a simple vertex-based mesh
        VertexMesh<2,2> mesh(4,6); // columns then rows

        // Set up cells
        std::vector<TissueCell> cells = SetUpCells(mesh);

        // Create tissue
        VertexBasedTissue<2> tissue(mesh, cells);

        // For coverage of WriteResultsToFiles()
        tissue.rGetCellUsingLocationIndex(0).SetCellType(TRANSIT);
        tissue.rGetCellUsingLocationIndex(0).SetMutationState(LABELLED);
        tissue.rGetCellUsingLocationIndex(1).SetCellType(DIFFERENTIATED);
        tissue.rGetCellUsingLocationIndex(1).SetMutationState(APC_ONE_HIT);
        tissue.rGetCellUsingLocationIndex(2).SetMutationState(APC_TWO_HIT);
        tissue.rGetCellUsingLocationIndex(3).SetMutationState(BETA_CATENIN_ONE_HIT);
        tissue.rGetCellUsingLocationIndex(4).SetCellType(APOPTOTIC);
        tissue.rGetCellUsingLocationIndex(4).StartApoptosis();
        tissue.rGetCellUsingLocationIndex(5).SetCellType(STEM);
        tissue.SetCellAncestorsToNodeIndices();

        std::string output_directory = "TestVertexBasedTissueOutputWriters";
        OutputFileHandler output_file_handler(output_directory, false);

        TS_ASSERT_THROWS_NOTHING(tissue.CreateOutputFiles(output_directory, false, true, true, false, true, true));

        tissue.WriteResultsToFiles(true, true, false, true, true);

        TS_ASSERT_THROWS_NOTHING(tissue.CloseOutputFiles(true, true, false, true, true));

        // Compare output with saved files of what they should look like
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.viznodes     notforrelease-cancer/test/data/TestVertexBasedTissueOutputWriters/results.viznodes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizelements     notforrelease-cancer/test/data/TestVertexBasedTissueOutputWriters/results.vizelements").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizcelltypes     notforrelease-cancer/test/data/TestVertexBasedTissueOutputWriters/results.vizcelltypes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizancestors     notforrelease-cancer/test/data/TestVertexBasedTissueOutputWriters/results.vizancestors").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellmutationstates.dat     notforrelease-cancer/test/data/TestVertexBasedTissueOutputWriters/cellmutationstates.dat").c_str()), 0);

        // For coverage
        TS_ASSERT_THROWS_NOTHING(tissue.WriteResultsToFiles(true, false, false, true, false));
    }

    // At the moment the tissue cannot be properly archived since the mesh cannot be. This test
    // just checks that the cells are correctly archived.
    void TestArchivingVertexBasedTissue() throw (Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "tissue.arch";

        // Archive tissue
        {
            // Need to set up time
            unsigned num_steps=10;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);

            // Create mesh
            VertexMeshReader2d mesh_reader("notforrelease-cancer/test/data/TestVertexMesh/vertex_mesh");
            VertexMesh<2,2> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            // Set up cells
            std::vector<TissueCell> cells = SetUpCells(mesh);

            // Create tissue
            VertexBasedTissue<2>* const p_tissue = new VertexBasedTissue<2>(mesh, cells);

            // Cells have been given birth times of 0 and -1.
            // Loop over them to run to time 0.0;
            for (AbstractTissue<2>::Iterator cell_iter=p_tissue->Begin();
                 cell_iter!=p_tissue->End();
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

            // Tidy up
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

            // Restore the tissue
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> *p_simulation_time;

            VertexBasedTissue<2>* p_tissue;

            // The following line is required because the loading of a tissue
            // is usually called by the method TissueSimulation::Load()
            MeshArchiveInfo::meshPathname = "notforrelease-cancer/test/data/TestVertexMesh/vertex_mesh";

            input_arch >> p_tissue;

            // Cells have been given birth times of 0, -1, -2, -3, -4.
            // this checks that individual cells and their models are archived.
            unsigned counter = 0u;
            for (AbstractTissue<2>::Iterator cell_iter=p_tissue->Begin();
                 cell_iter!=p_tissue->End();
                 ++cell_iter)
            {
                TS_ASSERT_DELTA(cell_iter->GetAge(), (double)(counter), 1e-7);
                counter++;
            }

            // Check the simulation time has been restored (through the cell)
            TS_ASSERT_EQUALS(p_simulation_time->GetTime(), 0.0);

            // Check the tissue has been restored
            TS_ASSERT_EQUALS(p_tissue->rGetCells().size(), 2u);

            // Tidy up
            delete p_tissue;
        }
    }
    
    void TestUpdateNodeLocations()
    {
        // Create a simple 2D VertexMesh
        VertexMesh<2,2> mesh(5,3); // columns then rows

        // Set up cells
        std::vector<TissueCell> cells = SetUpCells(mesh);
        
        // Create tissue
        VertexBasedTissue<2> tissue(mesh, cells);
        
        // Make up some forces
        std::vector<c_vector<double, 2> > old_posns(tissue.GetNumNodes());
        std::vector<c_vector<double, 2> > forces_on_nodes(tissue.GetNumNodes());
                
        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
            old_posns[i][0] = mesh.GetNode(i)->rGetLocation()[0];
            old_posns[i][1] = mesh.GetNode(i)->rGetLocation()[1];

            forces_on_nodes[i][0] = i*0.01;
            forces_on_nodes[i][1] = 2*i*0.01;
        }
        
        double time_step = 0.01;

        tissue.UpdateNodeLocations(forces_on_nodes, time_step);
       
        for (unsigned i=0; i<tissue.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA(tissue.GetNode(i)->rGetLocation()[0], old_posns[i][0] +   i*0.01*0.01, 1e-9);
            TS_ASSERT_DELTA(tissue.GetNode(i)->rGetLocation()[1], old_posns[i][1] + 2*i*0.01*0.01, 1e-9);
        }
    }
};


#endif /*TESTVERTEXBASEDTISSUE_HPP_*/
