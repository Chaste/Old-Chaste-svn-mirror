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
#ifndef TESTBOXCOLLECTION_HPP_
#define TESTBOXCOLLECTION_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "TetrahedralMesh.hpp"
#include "BoxCollection.hpp"
#include "TrianglesMeshReader.hpp"

class TestBoxCollection : public CxxTest::TestSuite
{
public:
    void TestBox() throw (Exception)
    {
        c_vector<double, 2*2> box_size;
        box_size(0) = -0.1; // min x
        box_size(1) = 1.1; // max x
        box_size(2) = -0.1; // min y
        box_size(3) = 1.1; // max y

        Box<2> test_box(box_size);
        c_vector<double, 2*2> returned_min_max_values = test_box.rGetMinAndMaxValues();
        for (unsigned i=0; i<4; i++)
        {
            TS_ASSERT_EQUALS(returned_min_max_values(i), box_size(i));
        }

        c_vector<double, 2> node_location;
        node_location(0) = 0.5;
        node_location(1) = 0.5;

        Node<2> test_node(213, node_location);

        test_box.AddNode(&test_node);
        std::set< Node<2>* > nodes_contained_before = test_box.rGetNodesContained();

        TS_ASSERT_EQUALS(*(nodes_contained_before.begin()), &test_node);
        TS_ASSERT_EQUALS((*(nodes_contained_before.begin()))->GetIndex(), 213u);

        test_box.RemoveNode(&test_node);
        std::set< Node<2>* > nodes_contained_after = test_box.rGetNodesContained();
        TS_ASSERT(nodes_contained_after.empty());
    }


    void TestBoxGeneration1d() throw (Exception)
    {
        // Create a mesh
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructLinearMesh(20);

        double cut_off_length = 5.0;
        c_vector<double, 2> domain_size;
        domain_size(0) = -0.1;
        domain_size(1) = 20.15;

        BoxCollection<1> box_collection(cut_off_length, domain_size);
        
        box_collection.SetupAllLocalBoxes();

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            unsigned box_index = box_collection.CalculateContainingBox(mesh.GetNode(i));
            box_collection.rGetBox(box_index).AddNode(mesh.GetNode(i));
        }

        TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), 5u);

        for (unsigned i=0; i<box_collection.GetNumBoxes(); i++)
        {
            std::set< Node<1>* > nodes_in_box = box_collection.rGetBox(i).rGetNodesContained();
            c_vector<double, 2> box_min_max_values = box_collection.rGetBox(i).rGetMinAndMaxValues();

            for (std::set< Node<1>* >::iterator it_nodes_in_box = nodes_in_box.begin();
                 it_nodes_in_box != nodes_in_box.end();
                 it_nodes_in_box++)
            {
                Node<1>* current_node = *it_nodes_in_box;
                double x_position = current_node->rGetLocation()[0];

                double epsilon = 1e-12;

                TS_ASSERT_LESS_THAN(box_min_max_values(0)-epsilon, x_position);
                TS_ASSERT_LESS_THAN(x_position, box_min_max_values(1)+epsilon);
            }
        }

        std::set<unsigned> local_boxes_to_box_0 = box_collection.GetLocalBoxes(0);
        std::set<unsigned> correct_answer_0;
        correct_answer_0.insert(0);
        correct_answer_0.insert(1);
        TS_ASSERT_EQUALS(local_boxes_to_box_0, correct_answer_0);

        std::set<unsigned> local_boxes_to_box_1 = box_collection.GetLocalBoxes(1);
        std::set<unsigned> correct_answer_1;
        correct_answer_1.insert(0);
        correct_answer_1.insert(1);
        correct_answer_1.insert(2);
        TS_ASSERT_EQUALS(local_boxes_to_box_1, correct_answer_1);

        std::set<unsigned> local_boxes_to_box_4 = box_collection.GetLocalBoxes(4);
        std::set<unsigned> correct_answer_4;
        correct_answer_4.insert(3);
        correct_answer_4.insert(4);
        TS_ASSERT_EQUALS(local_boxes_to_box_4, correct_answer_4);

        c_vector<double,1> miles_away;
        miles_away(0) = 47323854;
        TS_ASSERT_THROWS_CONTAINS(box_collection.CalculateContainingBox(miles_away), "The point provided in outside all of the boxes");
    }

    // very simple test
    void TestAddElement() throw(Exception)
    {
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructRegularSlabMesh(0.5, 1.0);
 
        double width = 0.4;
        c_vector<double, 2> domain_size;
        domain_size(0) = 0.0;
        domain_size(1) = 1.0;

        BoxCollection<1> box_collection(width, domain_size);
        box_collection.rGetBox(0).AddElement(mesh.GetElement(0));
        TS_ASSERT_EQUALS(box_collection.rGetBox(0).rGetElementsContained().size(), 1u);
        TS_ASSERT_EQUALS(box_collection.rGetBox(1).rGetElementsContained().size(), 0u);
        TS_ASSERT_EQUALS(box_collection.rGetBox(2).rGetElementsContained().size(), 0u);
        TS_ASSERT_EQUALS(*(box_collection.rGetBox(0).rGetElementsContained().begin()), mesh.GetElement(0));
    }
    

    void TestSetupAllLocalBoxes2d() throw(Exception)
    {
        double width = 1.0;

        c_vector<double, 2*2> domain_size;
        domain_size(0) = 0;
        domain_size(1) = 4-0.01;
        domain_size(2) = 0;
        domain_size(3) = 3-0.01;

        BoxCollection<2> box_collection(width, domain_size);

        assert(box_collection.GetNumBoxes()==12); // 4 * 3 boxes altogether

        box_collection.SetupAllLocalBoxes();
        
        std::set<unsigned> local_boxes_to_box_0 = box_collection.GetLocalBoxes(0);

        std::set<unsigned> correct_answer_0;
        correct_answer_0.insert(0);
        correct_answer_0.insert(1);
        correct_answer_0.insert(4);
        correct_answer_0.insert(5);
        TS_ASSERT_EQUALS(local_boxes_to_box_0, correct_answer_0);

        std::set<unsigned> local_boxes_to_box_3 = box_collection.GetLocalBoxes(3);
        std::set<unsigned> correct_answer_3;
        correct_answer_3.insert(3);
        correct_answer_3.insert(2);
        correct_answer_3.insert(6);
        correct_answer_3.insert(7);
        TS_ASSERT_EQUALS(local_boxes_to_box_3, correct_answer_3);

        std::set<unsigned> local_boxes_to_box_5 = box_collection.GetLocalBoxes(5);
        std::set<unsigned> correct_answer_5;
        correct_answer_5.insert(0);
        correct_answer_5.insert(1);
        correct_answer_5.insert(2);
        correct_answer_5.insert(4);
        correct_answer_5.insert(5);
        correct_answer_5.insert(6);
        correct_answer_5.insert(8);
        correct_answer_5.insert(9);
        correct_answer_5.insert(10);
        TS_ASSERT_EQUALS(local_boxes_to_box_5, correct_answer_5);

        std::set<unsigned> local_boxes_to_box_10 = box_collection.GetLocalBoxes(10);
        std::set<unsigned> correct_answer_10;
        correct_answer_10.insert(5);
        correct_answer_10.insert(6);
        correct_answer_10.insert(7);
        correct_answer_10.insert(9);
        correct_answer_10.insert(10);
        correct_answer_10.insert(11);
        TS_ASSERT_EQUALS(local_boxes_to_box_10, correct_answer_10);
    }



    void TestSetupAllLocalBoxes3d() throw(Exception)
    {
        double width = 1.0;

        c_vector<double, 2*3> domain_size;
        domain_size(0) = 0;
        domain_size(1) = 4-0.01;
        domain_size(2) = 0;
        domain_size(3) = 3-0.01;
        domain_size(4) = 0;
        domain_size(5) = 2-0.01;

        BoxCollection<3> box_collection(width, domain_size);

        assert(box_collection.GetNumBoxes()==24); // 4 * 3 * 2 boxes altogether

        box_collection.SetupAllLocalBoxes();
        
        std::set<unsigned> local_boxes_to_box_0 = box_collection.GetLocalBoxes(0);

        std::set<unsigned> correct_answer_0;
        correct_answer_0.insert(0);
        correct_answer_0.insert(1);
        correct_answer_0.insert(4);
        correct_answer_0.insert(5);
        correct_answer_0.insert(12);
        correct_answer_0.insert(13);
        correct_answer_0.insert(16);
        correct_answer_0.insert(17);
        TS_ASSERT_EQUALS(local_boxes_to_box_0, correct_answer_0);

        std::set<unsigned> local_boxes_to_box_3 = box_collection.GetLocalBoxes(3);
        std::set<unsigned> correct_answer_3;
        correct_answer_3.insert(3);
        correct_answer_3.insert(2);
        correct_answer_3.insert(6);
        correct_answer_3.insert(7);
        correct_answer_3.insert(14);
        correct_answer_3.insert(15);
        correct_answer_3.insert(18);
        correct_answer_3.insert(19);
        TS_ASSERT_EQUALS(local_boxes_to_box_3, correct_answer_3);

        std::set<unsigned> local_boxes_to_box_5 = box_collection.GetLocalBoxes(5);
        std::set<unsigned> correct_answer_5;
        correct_answer_5.insert(0);
        correct_answer_5.insert(1);
        correct_answer_5.insert(2);
        correct_answer_5.insert(4);
        correct_answer_5.insert(5);
        correct_answer_5.insert(6);
        correct_answer_5.insert(8);
        correct_answer_5.insert(9);
        correct_answer_5.insert(10);
        correct_answer_5.insert(12);
        correct_answer_5.insert(13);
        correct_answer_5.insert(14);
        correct_answer_5.insert(16);
        correct_answer_5.insert(17);
        correct_answer_5.insert(18);
        correct_answer_5.insert(20);
        correct_answer_5.insert(21);
        correct_answer_5.insert(22);

        TS_ASSERT_EQUALS(local_boxes_to_box_5, correct_answer_5);

        std::set<unsigned> local_boxes_to_box_19 = box_collection.GetLocalBoxes(19);
        std::set<unsigned> correct_answer_19;
        correct_answer_19.insert(2);
        correct_answer_19.insert(3);
        correct_answer_19.insert(6);
        correct_answer_19.insert(7);
        correct_answer_19.insert(10);
        correct_answer_19.insert(11);
        correct_answer_19.insert(14);
        correct_answer_19.insert(15);
        correct_answer_19.insert(18);
        correct_answer_19.insert(19);
        correct_answer_19.insert(22);
        correct_answer_19.insert(23);
        TS_ASSERT_EQUALS(local_boxes_to_box_19, correct_answer_19);
        
        std::set<unsigned> local_boxes_to_box_22 = box_collection.GetLocalBoxes(22);
        std::set<unsigned> correct_answer_22;
        correct_answer_22.insert(5);
        correct_answer_22.insert(6);
        correct_answer_22.insert(7);
        correct_answer_22.insert(9);
        correct_answer_22.insert(10);
        correct_answer_22.insert(11);
        correct_answer_22.insert(17);
        correct_answer_22.insert(18);
        correct_answer_22.insert(19);
        correct_answer_22.insert(21);
        correct_answer_22.insert(22);
        correct_answer_22.insert(23);
        TS_ASSERT_EQUALS(local_boxes_to_box_22, correct_answer_22);
        
    }

    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
    //
    //  Cancer/cell_based tests
    //  The following are tests written from this BoxCollection used to be
    //  cell_based/src/tissue/NodeBoxCollection and test the cancer related
    //  functionality
    //
    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
    void TestPairsReturned1d() throw (Exception)
    {
        std::vector< ChastePoint<1>* > points(5);
        points[0] = new ChastePoint<1>(0.2);
        points[1] = new ChastePoint<1>(0.7);
        points[2] = new ChastePoint<1>(1.3);
        points[3] = new ChastePoint<1>(2.1);
        points[4] = new ChastePoint<1>(6.9);

        std::vector<Node<1>* > nodes;
        for (unsigned i=0; i<points.size(); i++)
        {
            nodes.push_back(new Node<1>(i, *(points[i]), false));
        }

        double cut_off_length = 1.0;

        c_vector<double, 2> domain_size;
        domain_size(0) = 0.0;
        domain_size(1) = 7.0;

        BoxCollection<1> box_collection(cut_off_length, domain_size);
        
        box_collection.SetupLocalBoxesHalfOnly();

        for (unsigned i=0; i<nodes.size(); i++)
        {
            unsigned box_index = box_collection.CalculateContainingBox(nodes[i]);
            box_collection.rGetBox(box_index).AddNode(nodes[i]);
        }

        std::set< std::pair<Node<1>*, Node<1>* > > pairs_returned;
        box_collection.CalculateNodePairs(nodes,pairs_returned);

        std::set< std::pair<Node<1>*, Node<1>* > > pairs_should_be;
        pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[0],nodes[1]));
        pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[0],nodes[2]));
        pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[1],nodes[2]));
        pairs_should_be.insert(std::pair<Node<1>*, Node<1>*>(nodes[2],nodes[3]));

        TS_ASSERT_EQUALS(pairs_should_be, pairs_returned);

        for (unsigned i=0; i<points.size(); i++)
        {
            delete nodes[i];
            delete points[i];
        }
    }

    void TestBoxGeneration2d() throw (Exception)
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        double cut_off_length = 0.2;

        c_vector<double, 2*2> domain_size;
        domain_size(0) = -0.1;
        domain_size(1) = 1.15;
        domain_size(2) = -0.1;
        domain_size(3) = 1.15;

        BoxCollection<2> box_collection(cut_off_length, domain_size);

        box_collection.SetupLocalBoxesHalfOnly();

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            unsigned box_index = box_collection.CalculateContainingBox(mesh.GetNode(i));
            box_collection.rGetBox(box_index).AddNode(mesh.GetNode(i));
        }

        TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), 49u);

        for (unsigned i=0; i<box_collection.GetNumBoxes(); i++)
        {
            std::set< Node<2>* > nodes_in_box = box_collection.rGetBox(i).rGetNodesContained();
            c_vector<double, 2*2> box_min_max_values = box_collection.rGetBox(i).rGetMinAndMaxValues();

            for (std::set< Node<2>* >::iterator it_nodes_in_box = nodes_in_box.begin();
                 it_nodes_in_box != nodes_in_box.end();
                 it_nodes_in_box++)
            {
                Node<2>* current_node = *it_nodes_in_box;
                double x_position = current_node->rGetLocation()[0];
                double y_position = current_node->rGetLocation()[1];

                double epsilon = 1e-12;

                TS_ASSERT_LESS_THAN(box_min_max_values(0)-epsilon, x_position);
                TS_ASSERT_LESS_THAN(x_position, box_min_max_values(1)+epsilon);
                TS_ASSERT_LESS_THAN(box_min_max_values(2)-epsilon, y_position);
                TS_ASSERT_LESS_THAN(y_position, box_min_max_values(3)+epsilon);
            }
        }

        // Have checked that all the local boxes are calculated correctly on a 5 by 6 grid - here we
        // hardcode a few checks on the 7 by 7 grid.
        std::set<unsigned> local_boxes_to_box_0 = box_collection.GetLocalBoxes(0);
        std::set<unsigned> correct_answer_0;
        correct_answer_0.insert(0);
        correct_answer_0.insert(1);
        correct_answer_0.insert(7);
        correct_answer_0.insert(8);
        TS_ASSERT_EQUALS(local_boxes_to_box_0, correct_answer_0);

        std::set<unsigned> local_boxes_to_box_4 = box_collection.GetLocalBoxes(4);
        std::set<unsigned> correct_answer_4;
        correct_answer_4.insert(4);
        correct_answer_4.insert(5);
        correct_answer_4.insert(10);
        correct_answer_4.insert(11);
        correct_answer_4.insert(12);
        TS_ASSERT_EQUALS(local_boxes_to_box_4, correct_answer_4);

        std::set<unsigned> local_boxes_to_box_10 = box_collection.GetLocalBoxes(10);
        std::set<unsigned> correct_answer_10;
        correct_answer_10.insert(10);
        correct_answer_10.insert(11);
        correct_answer_10.insert(16);
        correct_answer_10.insert(17);
        correct_answer_10.insert(18);
        TS_ASSERT_EQUALS(local_boxes_to_box_10, correct_answer_10);

        std::set<unsigned> local_boxes_to_box_48 = box_collection.GetLocalBoxes(48);
        std::set<unsigned> correct_answer_48;
        correct_answer_48.insert(48);
        TS_ASSERT_EQUALS(local_boxes_to_box_48, correct_answer_48);
    }

    void TestPairsReturned2d() throw (Exception)
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

        double cut_off_length = 1.0;

        c_vector<double, 2*2> domain_size;
        domain_size(0) = 0.0;
        domain_size(1) = 4.0;
        domain_size(2) = 0.0;
        domain_size(3) = 4.0;

        BoxCollection<2> box_collection(cut_off_length, domain_size);

        box_collection.SetupLocalBoxesHalfOnly();


        for (unsigned i=0; i<nodes.size(); i++)
        {
            unsigned box_index = box_collection.CalculateContainingBox(nodes[i]);
            box_collection.rGetBox(box_index).AddNode(nodes[i]);
        }

        std::set< std::pair<Node<2>*, Node<2>* > > pairs_returned;
        box_collection.CalculateNodePairs(nodes,pairs_returned);

        std::set< std::pair<Node<2>*, Node<2>* > > pairs_should_be;
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[0],nodes[1]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[2],nodes[7]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[2],nodes[8]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[4]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[5]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[7]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[8]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[3],nodes[2]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[4],nodes[5]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[4],nodes[7]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[4],nodes[8]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[4],nodes[2]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[5],nodes[6]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[5],nodes[7]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[5],nodes[8]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[5],nodes[2]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[6],nodes[7]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[6],nodes[8]));
        pairs_should_be.insert(std::pair<Node<2>*, Node<2>*>(nodes[7],nodes[8]));

        TS_ASSERT_EQUALS(pairs_should_be, pairs_returned);

        for (unsigned i=0; i<points.size(); i++)
        {
            delete nodes[i];
            delete points[i];
        }
    }


    void TestBoxGeneration3d() throw (Exception)
    {
        // Create a mesh
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(4,5,6);

        double cut_off_length = 2.0;

        c_vector<double, 2*3> domain_size;
        domain_size(0) = -0.1;
        domain_size(1) = 4.15;
        domain_size(2) = -0.1;
        domain_size(3) = 5.15;
        domain_size(4) = -0.1;
        domain_size(5) = 6.15;

        BoxCollection<3> box_collection(cut_off_length, domain_size);

        box_collection.SetupLocalBoxesHalfOnly();

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            unsigned box_index = box_collection.CalculateContainingBox(mesh.GetNode(i));
            box_collection.rGetBox(box_index).AddNode(mesh.GetNode(i));
        }

        TS_ASSERT_EQUALS(box_collection.GetNumBoxes(), 36u);

        for (unsigned i=0; i<box_collection.GetNumBoxes(); i++)
        {
            std::set< Node<3>* > nodes_in_box = box_collection.rGetBox(i).rGetNodesContained();
            c_vector<double, 2*3> box_min_max_values = box_collection.rGetBox(i).rGetMinAndMaxValues();

            for (std::set< Node<3>* >::iterator it_nodes_in_box = nodes_in_box.begin();
                 it_nodes_in_box != nodes_in_box.end();
                 it_nodes_in_box++)
            {
                Node<3>* current_node = *it_nodes_in_box;
                double x_position = current_node->rGetLocation()[0];
                double y_position = current_node->rGetLocation()[1];
                double z_position = current_node->rGetLocation()[2];

                double epsilon = 1e-12;

                TS_ASSERT_LESS_THAN(box_min_max_values(0)-epsilon, x_position);
                TS_ASSERT_LESS_THAN(x_position, box_min_max_values(1)+epsilon);
                TS_ASSERT_LESS_THAN(box_min_max_values(2)-epsilon, y_position);
                TS_ASSERT_LESS_THAN(y_position, box_min_max_values(3)+epsilon);
                TS_ASSERT_LESS_THAN(box_min_max_values(4)-epsilon, z_position);
                TS_ASSERT_LESS_THAN(z_position, box_min_max_values(5)+epsilon);
            }
        }

        std::set<unsigned> local_boxes_to_box_0 = box_collection.GetLocalBoxes(0);
        std::set<unsigned> correct_answer_0;
        correct_answer_0.insert(0);
        correct_answer_0.insert(1);
        correct_answer_0.insert(3);
        correct_answer_0.insert(4);
        correct_answer_0.insert(9);
        correct_answer_0.insert(10);
        correct_answer_0.insert(12);
        correct_answer_0.insert(13);
        TS_ASSERT_EQUALS(local_boxes_to_box_0, correct_answer_0);

        std::set<unsigned> local_boxes_to_box_13 = box_collection.GetLocalBoxes(13);
        std::set<unsigned> correct_answer_13;
        correct_answer_13.insert(4);
        correct_answer_13.insert(5);
        correct_answer_13.insert(6);
        correct_answer_13.insert(7);
        correct_answer_13.insert(8);
        correct_answer_13.insert(13);
        correct_answer_13.insert(14);
        correct_answer_13.insert(15);
        correct_answer_13.insert(16);
        correct_answer_13.insert(17);
        correct_answer_13.insert(22);
        correct_answer_13.insert(23);
        correct_answer_13.insert(24);
        correct_answer_13.insert(25);
        correct_answer_13.insert(26);
        TS_ASSERT_EQUALS(local_boxes_to_box_13, correct_answer_13);

        std::set<unsigned> local_boxes_to_box_34 = box_collection.GetLocalBoxes(34);
        std::set<unsigned> correct_answer_34;
        correct_answer_34.insert(25);
        correct_answer_34.insert(26);
        correct_answer_34.insert(34);
        correct_answer_34.insert(35);
        TS_ASSERT_EQUALS(local_boxes_to_box_34, correct_answer_34);

        std::set<unsigned> local_boxes_to_box_35 = box_collection.GetLocalBoxes(35);
        std::set<unsigned> correct_answer_35;
        correct_answer_35.insert(26);
        correct_answer_35.insert(35);
        TS_ASSERT_EQUALS(local_boxes_to_box_35, correct_answer_35);
    }
};

#endif /*TESTBOXCOLLECTION_HPP_*/
