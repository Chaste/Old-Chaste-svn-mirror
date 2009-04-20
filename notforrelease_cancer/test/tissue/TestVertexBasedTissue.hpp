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

#ifndef TESTVERTEXBASEDTISSUE_HPP_
#define TESTVERTEXBASEDTISSUE_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "VertexBasedTissue.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
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
            TissueCell cell(DIFFERENTIATED, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
            double birth_time = 0.0 - i;
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        return cells;
    }

public:

    // Test construction, accessors and iterator
    void TestCreateSmallVertexBasedTissue() throw (Exception)
    {
        // Create a simple 2D VertexMesh
        VertexMesh<2,2> mesh(5, 3, 0.01, 2.0);

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
            TS_ASSERT_EQUALS(tissue.GetLocationIndexUsingCell(&(*cell_iter)), counter);

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


    void TestValidate() throw (Exception)
    {
        // Create a simple vertex-based mesh
        VertexMesh<2,2> mesh(3, 3, 0.01, 2.0);

        // Set up cells, one for each element.
        // Get each a birth time of -element_index, so the age = element_index.
        std::vector<TissueCell> cells;
        std::vector<unsigned> cell_location_indices;
        for (unsigned i=0; i<mesh.GetNumElements()-1; i++)
        {
            TissueCell cell(STEM, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
            double birth_time = 0.0 - i;
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
            cell_location_indices.push_back(i);
        }

        // This should throw an exception as the number of cells
        // does not equal the number of elements
        TS_ASSERT_THROWS_ANYTHING(VertexBasedTissue<2> tissue(mesh, cells));

        TissueCell cell(STEM, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
        double birth_time = 0.0 - mesh.GetNumElements()-1;
        cell.SetBirthTime(birth_time);
        cells.push_back(cell);
        cell_location_indices.push_back(mesh.GetNumElements()-1);

        // This should pass as the number of cells equals the number of elements
        TS_ASSERT_THROWS_NOTHING(VertexBasedTissue<2> tissue(mesh, cells));

        // Create tissue
        VertexBasedTissue<2> tissue(mesh, cells);

        // Check correspondence between elements and cells

        for (unsigned elem_index=0; elem_index<tissue.GetNumElements(); elem_index++)
        {
            std::set<unsigned> expected_node_indices;
            VertexElement<2,2>* p_expected_element = tissue.GetElement(elem_index);
            unsigned expected_index = p_expected_element->GetIndex();

            for (unsigned i=0; i<p_expected_element->GetNumNodes(); i++)
            {
                expected_node_indices.insert(p_expected_element->GetNodeGlobalIndex(i));
            }

            std::set<unsigned> actual_node_indices;
            TissueCell& r_cell = tissue.rGetCellUsingLocationIndex(elem_index);
            VertexElement<2,2>* p_actual_element = tissue.GetElementCorrespondingToCell(&r_cell);
            unsigned actual_index = p_actual_element->GetIndex();

            for (unsigned i=0; i<p_actual_element->GetNumNodes(); i++)
            {
                actual_node_indices.insert(p_actual_element->GetNodeGlobalIndex(i));
            }

            TS_ASSERT_EQUALS(actual_index, expected_index);
            TS_ASSERT_EQUALS(actual_node_indices, expected_node_indices);
        }
    }


    void TestGetTargetAreaOfCell() throw (Exception)
    {
        double apoptosis_time = CancerParameters::Instance()->GetApoptosisTime();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(apoptosis_time, 2);
        
        // Create mesh
        VertexMesh<2,2> mesh(3, 3, 0.01, 2.0);

        // Set up cells
        std::vector<TissueCell> cells;
        std::vector<unsigned> cell_location_indices;
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            CellType cell_type = STEM;

            if ((i==0) || (i==4))
            {
                cell_type = DIFFERENTIATED;
            }                       

            TissueCell cell(cell_type, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
            double birth_time = 0.0 - 2*i;
           
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);

            cell_location_indices.push_back(i);
        }

        // Create tissue
        VertexBasedTissue<2> tissue(mesh, cells);
        tissue.InitialiseCells(); // this method must be called explicitly as there is no simulation

        // Check GetTargetAreaOfCell()
        for (unsigned elem_index=0; elem_index<tissue.GetNumElements(); elem_index++)
        {
            TissueCell cell = tissue.rGetCellUsingLocationIndex(elem_index);
            double expected_area = CancerParameters::Instance()->GetMatureCellTargetArea();

            if (elem_index!=4 && elem_index<=7u)
            {
                expected_area *= 0.5*(1 + ((double)elem_index)/7.0);
            }

            double actual_area = tissue.GetTargetAreaOfCell(cell);

            TS_ASSERT_DELTA(actual_area, expected_area, 1e-12);
        }
        
        // Make 1 and 4 undergo apoptosis
        
        TissueCell cell_0 = tissue.rGetCellUsingLocationIndex(0);        
        TissueCell cell_1 = tissue.rGetCellUsingLocationIndex(1);
        TissueCell cell_4 = tissue.rGetCellUsingLocationIndex(4);
        
        cell_1.StartApoptosis();
        cell_4.StartApoptosis();
        
        double actual_area_0 = tissue.GetTargetAreaOfCell(cell_0);
        double actual_area_1 = tissue.GetTargetAreaOfCell(cell_1);
        double actual_area_4 = tissue.GetTargetAreaOfCell(cell_4);
                         
        double expected_area_0 = 0.5;                
                      
        double expected_area_1 = CancerParameters::Instance()->GetMatureCellTargetArea();
        expected_area_1 *= 0.5*(1.0 + 1.0/7.0); 
        
        double expected_area_4 = CancerParameters::Instance()->GetMatureCellTargetArea();
        
        TS_ASSERT_DELTA(actual_area_0, expected_area_0, 1e-12);
        TS_ASSERT_DELTA(actual_area_1, expected_area_1, 1e-12);
        TS_ASSERT_DELTA(actual_area_4, expected_area_4, 1e-12);
        
        p_simulation_time->IncrementTimeOneStep();
        
        double actual_area_0_after_dt = tissue.GetTargetAreaOfCell(cell_0);
        double actual_area_1_after_dt = tissue.GetTargetAreaOfCell(cell_1);
        double actual_area_4_after_dt = tissue.GetTargetAreaOfCell(cell_4);
        
        // Have run on for half the apoptosis time therefore the target area should have halved
        
        expected_area_0 = CancerParameters::Instance()->GetMatureCellTargetArea();
        expected_area_0 *= 0.5*(1.0 + 0.5*apoptosis_time/2.0); 
        
        TS_ASSERT_DELTA(actual_area_0_after_dt, expected_area_0, 1e-12);
        TS_ASSERT_DELTA(actual_area_1_after_dt, 0.5*expected_area_1, 1e-12);
        TS_ASSERT_DELTA(actual_area_4_after_dt, 0.5*expected_area_4, 1e-12);
        
        cell_0.StartApoptosis();
                
        // Now running on for a further half, i.e. the entire apoptosis time
        
        p_simulation_time->IncrementTimeOneStep();
        
        double actual_area_0_after_2dt = tissue.GetTargetAreaOfCell(cell_0);
        double actual_area_1_after_2dt = tissue.GetTargetAreaOfCell(cell_1);
        double actual_area_4_after_2dt = tissue.GetTargetAreaOfCell(cell_4);
        
        // Have run on for the further half of the apoptosis time therefore the target area 
        // should have gone to zero
        
        TS_ASSERT_DELTA(actual_area_0_after_2dt, 0.5*expected_area_0, 1e-12);
        TS_ASSERT_DELTA(actual_area_1_after_2dt, 0.0, 1e-12);
        TS_ASSERT_DELTA(actual_area_4_after_2dt, 0.0, 1e-12);                                    
               
    }


    void TestDampingConstant()
    {
        // Create mesh
        VertexMesh<2,2> mesh(3, 3, 0.01, 2.0);

        // Set up cells
        std::vector<TissueCell> cells;
        std::vector<unsigned> cell_location_indices;
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            TissueCell cell(STEM, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
            double birth_time = 0.0 - i;
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
            cell_location_indices.push_back(i);
        }
        cells[0].SetMutationState(APC_TWO_HIT);

        // Create tissue
        VertexBasedTissue<2> tissue(mesh, cells);
        tissue.InitialiseCells(); // this method must be called explicitly as there is no simulation

        // Test GetDampingConstant()
        double normal_damping_constant = CancerParameters::Instance()->GetDampingConstantNormal();
        double mutant_damping_constant = CancerParameters::Instance()->GetDampingConstantMutant();

        // Node 3 is contained in cell 2 only, therefore should have normal damping constant
        double damping_constant_at_node_3 = tissue.GetDampingConstant(3);
        TS_ASSERT_DELTA(damping_constant_at_node_3, normal_damping_constant, 1e-6);

        // Node 0 is contained in cell 0 only, therefore should have mutant damping constant
        double damping_constant_at_node_0 = tissue.GetDampingConstant(0);
        TS_ASSERT_DELTA(damping_constant_at_node_0, mutant_damping_constant, 1e-6);

        // Node 5 is contained in cells 0 and 1, therefore should an averaged damping constant
        double damping_constant_at_node_5 = tissue.GetDampingConstant(5);
        TS_ASSERT_DELTA(damping_constant_at_node_5, (normal_damping_constant+mutant_damping_constant)/2.0, 1e-6);

        // Node 9 is contained in cells 0, 1, 3, therefore should an averaged damping constant
        double damping_constant_at_node_9 = tissue.GetDampingConstant(9);
        TS_ASSERT_DELTA(damping_constant_at_node_9, (2*normal_damping_constant+mutant_damping_constant)/3.0, 1e-6);
    }


    void TestUpdateWithoutBirthOrDeath() throw (Exception)
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);

        // Create a simple vertex-based mesh
        VertexMesh<2,2> mesh(4, 6, 0.01, 2.0);

        // Set up cells
        std::vector<TissueCell> cells = SetUpCells(mesh);

        // Create tissue
        VertexBasedTissue<2> tissue(mesh, cells);

        unsigned num_cells_removed = tissue.RemoveDeadCells();
        TS_ASSERT_EQUALS(num_cells_removed, 0u);

        p_simulation_time->IncrementTimeOneStep();

        TS_ASSERT_THROWS_NOTHING(tissue.Update());
    }


    void TestAddCellWithSimpleMesh() throw (Exception)
    {
        // Make some nodes
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

        // For coverage, test GetLocationOfCellCentre()

        // Cell 0 is a rectangle with centre of mass (0,0)
        c_vector<double, 2> cell0_location = tissue.GetLocationOfCellCentre(&(tissue.rGetCellUsingLocationIndex(0)));
        TS_ASSERT_DELTA(cell0_location[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell0_location[1], 0.0, 1e-4);

        // Cell 1 is a triangle with centre of mass (0,4/3)
        c_vector<double, 2> cell1_location = tissue.GetLocationOfCellCentre(&(tissue.rGetCellUsingLocationIndex(1)));
        TS_ASSERT_DELTA(cell1_location[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell1_location[1], 4.0/3.0, 1e-4);

        unsigned old_num_nodes = vertex_mesh.GetNumNodes();
        unsigned old_num_elements = vertex_mesh.GetNumElements();
        unsigned old_num_cells = tissue.rGetCells().size();

        // Add new cell by dividing element 0 along short axis

        c_vector<double,2> new_cell_location = zero_vector<double>(2);

        TissueCell cell0 = tissue.rGetCellUsingLocationIndex(0);

        TissueCell new_cell(STEM, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
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
        TS_ASSERT_EQUALS(tissue.GetLocationIndexUsingCell(p_new_cell), old_num_elements);
    }


    void TestAddCellWithHoneycombMesh() throw (Exception)
    {
        // Create a mesh with 9 elements
        VertexMesh<2,2> vertex_mesh(3, 3, 0.01, 2.0);

        // Set up cells, one for each VertexElement. Give each cell
        // a random birth time of -elem_index, so its age is elem_index
        std::vector<TissueCell> cells;
        for (unsigned elem_index=0; elem_index<vertex_mesh.GetNumElements(); elem_index++)
        {
            CellType cell_type = DIFFERENTIATED;
            double birth_time = 0.0 - elem_index;

            // Cell 4 should divide immediately
            if (elem_index==4)
            {
                cell_type = STEM;
                birth_time = -50.0;
            }

            TissueCell cell(cell_type, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Create tissue
        VertexBasedTissue<2> tissue(vertex_mesh, cells);

        // Initialise cells (usually called in tissue simulation constructor)
        tissue.InitialiseCells();

        unsigned old_num_nodes = vertex_mesh.GetNumNodes();
        unsigned old_num_elements = vertex_mesh.GetNumElements();
        unsigned old_num_cells = tissue.rGetCells().size();

        // Add a new cell by dividing cell 4

        tissue.rGetCellUsingLocationIndex(4).ReadyToDivide();
        TissueCell& cell4 = tissue.rGetCellUsingLocationIndex(4);

        TissueCell new_cell = tissue.rGetCellUsingLocationIndex(4).Divide();

        c_vector<double, 2> new_location = zero_vector<double>(2);

        // Add new cell to the tissue
        tissue.AddCell(new_cell, new_location, &cell4);

        // Check that the new cell was successfully added to the tissue
        TS_ASSERT_EQUALS(tissue.GetNumNodes(), old_num_nodes+2);
        TS_ASSERT_EQUALS(tissue.GetNumElements(), old_num_elements+1);
        TS_ASSERT_EQUALS(tissue.rGetCells().size(), old_num_cells+1);
        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), old_num_elements+1);

        // Check the location of the new nodes
        TS_ASSERT_DELTA(tissue.GetNode(old_num_nodes)->rGetLocation()[0], 1.8509, 1e-4);
        TS_ASSERT_DELTA(tissue.GetNode(old_num_nodes)->rGetLocation()[1], 1.7058, 1e-4);

        TS_ASSERT_DELTA(tissue.GetNode(old_num_nodes+1)->rGetLocation()[0], 1.0358, 1e-4);
        TS_ASSERT_DELTA(tissue.GetNode(old_num_nodes+1)->rGetLocation()[1], 2.2941, 1e-4);

        // Now test the nodes in each element
        for (unsigned i=0; i<tissue.GetNumElements(); i++)
        {
            if (i==4 || i==9)
            {
                // Elements 4 and 9 should each have one less node
                TS_ASSERT_EQUALS(tissue.GetElement(i)->GetNumNodes(), 5u);
            }
            else if (i==5 || i==6)
            {
                // Elements 5 and 6 should each have one extra node
                TS_ASSERT_EQUALS(tissue.GetElement(i)->GetNumNodes(), 7u);
            }
            else
            {
                TS_ASSERT_EQUALS(tissue.GetElement(i)->GetNumNodes(), 6u);
            }
        }

        // Check node ownership for a few elements

        TS_ASSERT_EQUALS(tissue.GetElement(0)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(tissue.GetElement(0)->GetNodeGlobalIndex(1), 1u);
        TS_ASSERT_EQUALS(tissue.GetElement(0)->GetNodeGlobalIndex(2), 5u);
        TS_ASSERT_EQUALS(tissue.GetElement(0)->GetNodeGlobalIndex(3), 9u);
        TS_ASSERT_EQUALS(tissue.GetElement(0)->GetNodeGlobalIndex(4), 8u);
        TS_ASSERT_EQUALS(tissue.GetElement(0)->GetNodeGlobalIndex(5), 4u);

        TS_ASSERT_EQUALS(tissue.GetElement(4)->GetNodeGlobalIndex(0), 30u);
        TS_ASSERT_EQUALS(tissue.GetElement(4)->GetNodeGlobalIndex(1), 18u);
        TS_ASSERT_EQUALS(tissue.GetElement(4)->GetNodeGlobalIndex(2), 22u);
        TS_ASSERT_EQUALS(tissue.GetElement(4)->GetNodeGlobalIndex(3), 21u);
        TS_ASSERT_EQUALS(tissue.GetElement(4)->GetNodeGlobalIndex(4), 31u);

        TS_ASSERT_EQUALS(tissue.GetElement(5)->GetNodeGlobalIndex(0), 10u);
        TS_ASSERT_EQUALS(tissue.GetElement(5)->GetNodeGlobalIndex(1), 11u);
        TS_ASSERT_EQUALS(tissue.GetElement(5)->GetNodeGlobalIndex(2), 15u);
        TS_ASSERT_EQUALS(tissue.GetElement(5)->GetNodeGlobalIndex(3), 19u);
        TS_ASSERT_EQUALS(tissue.GetElement(5)->GetNodeGlobalIndex(4), 18u);
        TS_ASSERT_EQUALS(tissue.GetElement(5)->GetNodeGlobalIndex(5), 30u);
        TS_ASSERT_EQUALS(tissue.GetElement(5)->GetNodeGlobalIndex(6), 14u);

        // Test element ownership for a few nodes

        std::set<unsigned> expected_elements_containing_node_5;
        expected_elements_containing_node_5.insert(0);
        expected_elements_containing_node_5.insert(1);

        TS_ASSERT_EQUALS(tissue.GetNode(5)->rGetContainingElementIndices(), expected_elements_containing_node_5);

        std::set<unsigned> expected_elements_containing_node_13;
        expected_elements_containing_node_13.insert(1);
        expected_elements_containing_node_13.insert(3);
        expected_elements_containing_node_13.insert(9);

        TS_ASSERT_EQUALS(tissue.GetNode(13)->rGetContainingElementIndices(), expected_elements_containing_node_13);

        std::set<unsigned> expected_elements_containing_node_30;
        expected_elements_containing_node_30.insert(5);
        expected_elements_containing_node_30.insert(4);
        expected_elements_containing_node_30.insert(9);

        TS_ASSERT_EQUALS(tissue.GetNode(30)->rGetContainingElementIndices(), expected_elements_containing_node_30);
    }


    void TestRemoveDeadCellsAndUpdate() throw (Exception)
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);

        // Create a simple vertex-based mesh
        VertexMesh<2,2> mesh(4, 6, 0.01, 2.0);

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

        TS_ASSERT_EQUALS(num_cells_removed, 1u);

        // We should now have one less real cell, since one cell has been
        // marked as dead, so is skipped by the tissue iterator
        TS_ASSERT_EQUALS(tissue.GetNumRealCells(), 23u);

        /// \todo Need some more tests here, on the new number of elements/nodes

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
            unsigned element_index = tissue.GetLocationIndexUsingCell(&(*cell_iter));
            element_indices.insert(element_index);
        }

        TS_ASSERT_EQUALS(element_indices, expected_elem_indices);
    }


    void TestVertexBasedTissueOutputWriters() throw (Exception)
    {
        // Create a simple vertex-based mesh
        VertexMesh<2,2> mesh(4, 6, 0.01, 2.0);

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

        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.viznodes     notforrelease_cancer/test/data/TestVertexBasedTissueOutputWriters/results.viznodes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizelements     notforrelease_cancer/test/data/TestVertexBasedTissueOutputWriters/results.vizelements").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizcelltypes     notforrelease_cancer/test/data/TestVertexBasedTissueOutputWriters/results.vizcelltypes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizancestors     notforrelease_cancer/test/data/TestVertexBasedTissueOutputWriters/results.vizancestors").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellmutationstates.dat     notforrelease_cancer/test/data/TestVertexBasedTissueOutputWriters/cellmutationstates.dat").c_str()), 0);

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
            unsigned num_steps = 10;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);

            // Create mesh
            VertexMeshReader2d mesh_reader("notforrelease_cancer/test/data/TestVertexMesh/vertex_mesh");
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
            MeshArchiveInfo::meshPathname = "notforrelease_cancer/test/data/TestVertexMesh/vertex_mesh";

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


    void TestUpdateNodeLocations() throw (Exception)
    {
        // Create a simple 2D VertexMesh
        VertexMesh<2,2> mesh(5, 3, 0.01, 2.0);

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


    /**
     * Test that post-#878, WntConcentration copes with a VertexBasedTissue.
     * \todo When vertex-based tissue code is added to cancer folder, move this
     *       test to TestWntConcentration.hpp
     */
    void TestWntConcentrationWithVertexBasedTissue() throw(Exception)
    {
        // Make some nodes
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

        // Create cells
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<vertex_mesh.GetNumElements(); i++)
        {
            TissueCell cell(DIFFERENTIATED, HEALTHY, new WntCellCycleModel());
            double birth_time = 0.0 - i;
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Create tissue
        VertexBasedTissue<2> tissue(vertex_mesh, cells);

        // Set the top of this tissue, for the purposes of computing the WntConcentration
        CancerParameters::Instance()->SetCryptLength(4.0);

        // Set up an instance of the WntConcentration singleton object
        WntConcentration* p_wnt = WntConcentration::Instance();
        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), false);

        // Check that the singleton can be set up
        p_wnt->SetType(LINEAR);
        p_wnt->SetTissue(tissue);

        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), true);

        // Check that the singleton can be destroyed then recreated
        WntConcentration::Destroy();
        WntConcentration::Instance()->SetType(NONE);
        WntConcentration::Instance()->SetTissue(tissue);
        TS_ASSERT_EQUALS(WntConcentration::Instance()->IsWntSetUp(), false); // not fully set up now it is a NONE type

        WntConcentration::Destroy();
        WntConcentration::Instance()->SetType(LINEAR);
        WntConcentration::Instance()->SetTissue(tissue);

        p_wnt = WntConcentration::Instance();
        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), true); // set up again

        double wnt_at_cell0 = p_wnt->GetWntLevel(&(tissue.rGetCellUsingLocationIndex(0)));
        double wnt_at_cell1 = p_wnt->GetWntLevel(&(tissue.rGetCellUsingLocationIndex(1)));

        // We have set the top of the tissue to be 4, so the WntConcentration should decrease linearly
        // up the tissue, from one at height 0 to zero at height 4.

        // Cell 0 has centre of mass (0,0)
        TS_ASSERT_DELTA(wnt_at_cell0, 1.0, 1e-4);

        // Cell 1 has centre of mass (0, 4/3)
        TS_ASSERT_DELTA(wnt_at_cell1, 2.0/3.0, 1e-4);


    }
};


#endif /*TESTVERTEXBASEDTISSUE_HPP_*/
