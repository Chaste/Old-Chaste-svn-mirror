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
#ifndef TESTNODEBOXCOLLECTION_HPP_
#define TESTNODEBOXCOLLECTION_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "NodeBasedTissue.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "TrianglesMeshReader.hpp"
#include "TetrahedralMesh.hpp"
#include "AbstractCancerTestSuite.hpp"


class TestNodeBoxCollection : public AbstractCancerTestSuite
{
public:
     
    void TestNodeBox() throw (Exception)
    {
        c_vector<double, 2*2> box_size;
        box_size(0) = -0.1; // min x
        box_size(1) = 1.1; // max x
        box_size(2) = -0.1; // min y
        box_size(3) = 1.1; // max y
        
        NodeBox<2> test_box(box_size);
        c_vector<double, 2*2> returned_min_max_values = test_box.rGetMinAndMaxValues();
        for (unsigned i=0;i<4;i++)
        {
            TS_ASSERT_EQUALS(returned_min_max_values(i),box_size(i));
        }
    
        c_vector<double, 2> node_location;
        node_location(0) = 0.5;
        node_location(1) = 0.5;
        
        Node<2> test_node(213,node_location);
        
        test_box.AddNode(&test_node);
        std::set< Node<2>* > nodes_contained_before = test_box.rGetNodesContained();
        
        TS_ASSERT_EQUALS(*(nodes_contained_before.begin()), &test_node);
        TS_ASSERT_EQUALS((*(nodes_contained_before.begin()))->GetIndex(), 213u);
        
        test_box.RemoveNode(&test_node);
        std::set< Node<2>* > nodes_contained_after = test_box.rGetNodesContained();
        TS_ASSERT(nodes_contained_after.empty());
    }
    
    
    void TestBoxGeneration2d() throw (Exception)
    {
        
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node. Get each a birth time of -node_index,
        // so the age = node_index
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            TissueCell cell(STEM, HEALTHY, new FixedDurationGenerationBasedCellCycleModel());
            double birth_time = 0.0 - i;
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        // Create a tissue
        NodeBasedTissue<2> node_based_tissue(mesh, cells);
        
        double cut_off_length = 0.2;
        
        c_vector<double, 2*2> domain_size;
        domain_size(0) = -0.1;
        domain_size(1) = 1.15;
        domain_size(2) = -0.1;
        domain_size(3) = 1.15;
            
        node_based_tissue.SplitUpIntoBoxes(cut_off_length, domain_size);
        
        TS_ASSERT_EQUALS(node_based_tissue.GetNodeBoxCollection()->GetNumBoxes(),49u);
        
        for (unsigned i=0;i<node_based_tissue.GetNodeBoxCollection()->GetNumBoxes();i++)
        {
            std::set< Node<2>* > nodes_in_box = node_based_tissue.GetNodeBoxCollection()->rGetBox(i).rGetNodesContained();
            c_vector<double, 2*2> box_min_max_values = node_based_tissue.GetNodeBoxCollection()->rGetBox(i).rGetMinAndMaxValues();
            
            for (std::set< Node<2>* >::iterator it_nodes_in_box = nodes_in_box.begin();
                it_nodes_in_box != nodes_in_box.end();
                it_nodes_in_box++)
            {
                Node<2>* current_node = *it_nodes_in_box;
                double x_position = current_node->rGetLocation()[0];
                double y_position = current_node->rGetLocation()[1];
                
                double epsilon = 1e-12;
                
                TS_ASSERT_LESS_THAN(box_min_max_values(0)-epsilon,x_position);
                TS_ASSERT_LESS_THAN(x_position,box_min_max_values(1)+epsilon);
                TS_ASSERT_LESS_THAN(box_min_max_values(2)-epsilon,y_position);
                TS_ASSERT_LESS_THAN(y_position,box_min_max_values(3)+epsilon);
            }
        }
        
        // have checked that all the local boxes are calculated correctly on a 5 by 6 grid - here we
        // hardcode a few checks on the 7 by 7 grid. 
        std::set<unsigned> local_boxes_to_box_0 = node_based_tissue.GetNodeBoxCollection()->GetLocalBoxes(0);
        std::set<unsigned> correct_answer_0;
        correct_answer_0.insert(0);
        correct_answer_0.insert(1);
        correct_answer_0.insert(7);
        correct_answer_0.insert(8);
        TS_ASSERT_EQUALS(local_boxes_to_box_0, correct_answer_0);

        std::set<unsigned> local_boxes_to_box_4 = node_based_tissue.GetNodeBoxCollection()->GetLocalBoxes(4);
        std::set<unsigned> correct_answer_4;
        correct_answer_4.insert(3);
        correct_answer_4.insert(4);
        correct_answer_4.insert(5);
        correct_answer_4.insert(10);
        correct_answer_4.insert(11);
        correct_answer_4.insert(12);
        TS_ASSERT_EQUALS(local_boxes_to_box_4, correct_answer_4);

        std::set<unsigned> local_boxes_to_box_10 = node_based_tissue.GetNodeBoxCollection()->GetLocalBoxes(10);
        std::set<unsigned> correct_answer_10; 
        correct_answer_10.insert(2);
        correct_answer_10.insert(3);
        correct_answer_10.insert(4);
        correct_answer_10.insert(9);
        correct_answer_10.insert(10);
        correct_answer_10.insert(11);
        correct_answer_10.insert(16);
        correct_answer_10.insert(17);
        correct_answer_10.insert(18);
        TS_ASSERT_EQUALS(local_boxes_to_box_10, correct_answer_10);


        std::set<unsigned> local_boxes_to_box_48 = node_based_tissue.GetNodeBoxCollection()->GetLocalBoxes(48);
        std::set<unsigned> correct_answer_48;
        correct_answer_48.insert(40);
        correct_answer_48.insert(41);
        correct_answer_48.insert(47);
        correct_answer_48.insert(48);
        TS_ASSERT_EQUALS(local_boxes_to_box_48, correct_answer_48);
    }
    
    void TestPairsReturned() throw (Exception)
    {
        std::vector< ChastePoint<2>* > points(10);
        points[0] = new ChastePoint<2>(0.2, 3.7);
        points[1] = new ChastePoint<2>(0.5, 3.2);
        points[2] = new ChastePoint<2>(1.1, 1.99);
        points[3] = new ChastePoint<2>(1.3, 0.8);
        points[4] = new ChastePoint<2>(1.3, 0.3);
        points[5] = new ChastePoint<2>(2.2, 0.6);
        points[6] = new ChastePoint<2>(3.5, 0.2);
        points[7] = new ChastePoint<2>(2.6, 1.4);
        points[8] = new ChastePoint<2>(2.4, 1.5);
        points[9] = new ChastePoint<2>(3.3, 3.6);
        
        std::vector<Node<2>* > nodes;
        for (unsigned i=0; i<points.size(); i++)
        {
            nodes.push_back(new Node<2>(i, *(points[i]), false));
        }

        NodeBasedTissue<2> node_based_tissue(nodes);
        
        double cut_off_length = 1.0;
        
        c_vector<double, 2*2> domain_size;
        domain_size(0) = 0.0;
        domain_size(1) = 4.0;
        domain_size(2) = 0.0;
        domain_size(3) = 4.0;
        
        node_based_tissue.SplitUpIntoBoxes(cut_off_length, domain_size);
        
        std::set< std::pair<Node<2>*, Node<2>* > > pairs_returned;
        node_based_tissue.GetNodeBoxCollection()->CalculateNodePairs(nodes,pairs_returned);
        
        std::set< std::pair<Node<2>*, Node<2>* > > pairs_should_be;
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[0],nodes[1]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[2],nodes[3]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[2],nodes[4]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[2],nodes[5]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[2],nodes[7]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[2],nodes[8]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[4]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[5]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[7]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[8]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[4],nodes[5]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[4],nodes[7]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[4],nodes[8]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[5],nodes[6]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[5],nodes[7]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[5],nodes[8]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[6],nodes[7]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[6],nodes[8]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[7],nodes[8]));
        
        TS_ASSERT_EQUALS(pairs_should_be,pairs_returned );
        
        for (unsigned i=0; i<points.size(); i++)
        {   // Tissue deletes the nodes
            delete points[i];
        }
    }
};

#endif /*TESTNODEBOXCOLLECTION_HPP_*/
