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

#ifndef TESTVERTEXELEMENT_HPP_
#define TESTVERTEXELEMENT_HPP_

#include <cxxtest/TestSuite.h>

#include "VertexElement.hpp"
#include "Element.hpp"
#include "Debug.hpp"

class TestVertexElement : public CxxTest::TestSuite
{
public:

    void Test1dVertexElementIn2d()
    {
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 1.0));

        VertexElement<1,2> element(0, nodes);

        // Test RegisterWithNodes()
        element.RegisterWithNodes();

        for (unsigned node_index=0; node_index<element.GetNumNodes(); node_index++)
        {
            TS_ASSERT_EQUALS(element.GetNode(node_index)->GetNumContainingElements(), 1u);
        }

        TS_ASSERT_EQUALS(element.GetNumNodes(), 2u);
        TS_ASSERT_EQUALS(element.GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(element.GetNode(1)->GetIndex(), 1u);

        // Test GetNodeLocalIndex()
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(0), 0u);
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(1), 1u);

        // Test Face methods
        TS_ASSERT_EQUALS(element.GetNumFaces(), 0u);
        VertexElement<0,2>* p_face = element.GetFace(0);
        TS_ASSERT(!p_face);
        TS_ASSERT_EQUALS(element.FaceIsOrientatedClockwise(0), false);

        // Test UpdateNode()
        Node<2>* p_node_2 = new Node<2>(2, false, 1.2, 1.3);
        element.UpdateNode(0, p_node_2);

        TS_ASSERT_DELTA(element.GetNode(0)->rGetLocation()[0], 1.2, 1e-12);
        TS_ASSERT_DELTA(element.GetNode(0)->rGetLocation()[1], 1.3, 1e-12);

        // Test ResetIndex()
        TS_ASSERT_EQUALS(element.GetIndex(), 0u);
        element.ResetIndex(5);
        TS_ASSERT_EQUALS(element.GetIndex(), 5u);

        // Test DeleteNode() and AddNode()
        TS_ASSERT_EQUALS(element.GetNumNodes(), 2u);
        element.DeleteNode(1);
        TS_ASSERT_EQUALS(element.GetNumNodes(), 1u);

        Node<2>* p_node_3 = new Node<2>(3, false, 0.1, 0.4);
        element.AddNode(0, p_node_3);
        TS_ASSERT_EQUALS(element.GetNumNodes(), 2u);

        TS_ASSERT_DELTA(element.GetNode(0)->rGetLocation()[0], 1.2, 1e-12);
        TS_ASSERT_DELTA(element.GetNode(0)->rGetLocation()[1], 1.3, 1e-12);
        TS_ASSERT_DELTA(element.GetNode(1)->rGetLocation()[0], 0.1, 1e-12);
        TS_ASSERT_DELTA(element.GetNode(1)->rGetLocation()[1], 0.4, 1e-12);

        // Test GetNodeLocalIndex()
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(0), UINT_MAX);
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(1), UINT_MAX);

        // Test MarkAsDeleted()
        element.MarkAsDeleted();

        for (unsigned node_index=0; node_index<element.GetNumNodes(); node_index++)
        {
            TS_ASSERT_EQUALS(element.GetNode(node_index)->GetNumContainingElements(), 0u);
        }

        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
        delete p_node_2;
        delete p_node_3;
    }

    void TestCreateVertexElement()
    {
        // Make 8 nodes to assign to a cube element
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, false, 0.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(4, false, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(5, false, 0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(6, false, 1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(7, false, 1.0, 1.0, 1.0));

        std::vector<Node<3>*> nodes_face_0, nodes_face_1, nodes_face_2, nodes_face_3, nodes_face_4, nodes_face_5;

        // Make 6 square faces out of these nodes
        nodes_face_0.push_back(nodes[0]);
        nodes_face_0.push_back(nodes[2]);
        nodes_face_0.push_back(nodes[4]);
        nodes_face_0.push_back(nodes[1]);

        nodes_face_1.push_back(nodes[4]);
        nodes_face_1.push_back(nodes[7]);
        nodes_face_1.push_back(nodes[5]);
        nodes_face_1.push_back(nodes[2]);

        nodes_face_2.push_back(nodes[7]);
        nodes_face_2.push_back(nodes[6]);
        nodes_face_2.push_back(nodes[1]);
        nodes_face_2.push_back(nodes[4]);

        nodes_face_3.push_back(nodes[0]);
        nodes_face_3.push_back(nodes[3]);
        nodes_face_3.push_back(nodes[5]);
        nodes_face_3.push_back(nodes[2]);

        nodes_face_4.push_back(nodes[1]);
        nodes_face_4.push_back(nodes[6]);
        nodes_face_4.push_back(nodes[3]);
        nodes_face_4.push_back(nodes[0]);

        nodes_face_5.push_back(nodes[7]);
        nodes_face_5.push_back(nodes[6]);
        nodes_face_5.push_back(nodes[3]);
        nodes_face_5.push_back(nodes[5]);

        std::vector<VertexElement<2,3>*> faces;
        faces.push_back(new VertexElement<2,3>(0, nodes_face_0));
        faces.push_back(new VertexElement<2,3>(1, nodes_face_1));
        faces.push_back(new VertexElement<2,3>(2, nodes_face_2));
        faces.push_back(new VertexElement<2,3>(3, nodes_face_3));
        faces.push_back(new VertexElement<2,3>(4, nodes_face_4));
        faces.push_back(new VertexElement<2,3>(5, nodes_face_5));

        std::vector<bool> orientations(faces.size());
        for (unsigned i=0; i<faces.size(); i++)
        {
            orientations[i] = true;
        }

        ///\todo Temporary test with hard-coded class

        // Make a cube element out of these faces
        VertexElement<3,3> element(0, faces, orientations);

        TS_ASSERT_EQUALS(element.GetNumNodes(),8u);
        TS_ASSERT_EQUALS(element.GetNumFaces(),6u);

        TS_ASSERT_EQUALS(element.GetIndex(),0u);

        // Test the position of some random nodes
        TS_ASSERT_DELTA(element.GetFace(0)->GetNode(0)->rGetLocation()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(element.GetFace(0)->GetNode(0)->rGetLocation()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(element.GetFace(0)->GetNode(0)->rGetLocation()[2], 0.0, 1e-6);

        TS_ASSERT_DELTA(element.GetFace(5)->GetNode(2)->rGetLocation()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(element.GetFace(5)->GetNode(2)->rGetLocation()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(element.GetFace(5)->GetNode(2)->rGetLocation()[2], 1.0, 1e-6);

        // Test orientations
        for (unsigned face_index=0; face_index<element.GetNumFaces(); face_index++)
        {
            TS_ASSERT_EQUALS(element.FaceIsOrientatedClockwise(face_index), true);
        }

        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
        for (unsigned i=0; i<faces.size(); i++)
        {
            delete faces[i];
        }
    }

    void TestVertexElementFaceConstructor()
    {
        // Create a regular hexagon
        std::vector<Node<2>*> nodes;
        std::vector<VertexElement<1,2>*> faces;
        std::vector<Node<2>*> face_nodes;
        std::vector<bool> orientations;
        unsigned num_nodes = 6;
        for (unsigned i=0; i<num_nodes; i++)
        {
            double theta = 2.0*M_PI*(double)(i)/(double)(num_nodes);
            nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
        }

        for (unsigned i=0; i<num_nodes-1; i++)
        {
            face_nodes.clear();
            face_nodes.push_back(nodes[i]);
            face_nodes.push_back(nodes[i+1]);

            faces.push_back(new VertexElement<1,2>(i, face_nodes));
            orientations.push_back(true);
        }

        // Create a face with negative orientation
        face_nodes.clear();
        face_nodes.push_back(nodes[0]);
        face_nodes.push_back(nodes[num_nodes-1]);
        faces.push_back(new VertexElement<1,2>(num_nodes-1, face_nodes));
        orientations.push_back(false);

        // Create element
        VertexElement<2,2> vertex_element(0, faces, orientations);

        TS_ASSERT_EQUALS(vertex_element.GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(vertex_element.GetNumFaces(), 6u);

        // Test that each face has the correct orientation and number of nodes
        for (unsigned face_index=0; face_index<6; face_index++)
        {
            TS_ASSERT_EQUALS(vertex_element.GetFace(face_index)->GetNumNodes(), 2u);

            bool is_clockwise = (face_index==5) ? false : true;
            TS_ASSERT_EQUALS(vertex_element.FaceIsOrientatedClockwise(face_index), is_clockwise);
        }

        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
            delete faces[i];
        }
    }

    void TestVertexElementDeleteAndAddNode()
    {
        // Create nodes and Faces
        std::vector<Node<2>*> nodes;
        unsigned num_nodes = 6;
        for (unsigned i=0; i<num_nodes; i++)
        {
            double theta = 2.0*M_PI*(double)(i)/(double)(num_nodes);
            nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
        }

        // Create element
        VertexElement<2,2> vertex_element(0, nodes);

        TS_ASSERT_EQUALS(vertex_element.GetNumNodes(), 6u);

        vertex_element.DeleteNode(3); // Removes (-1,0) node
        vertex_element.DeleteNode(0); // Removes (1,0) node

        // Test node is removed
        TS_ASSERT_EQUALS(vertex_element.GetNumNodes(), 4u);

        // Test other nodes are updated
        TS_ASSERT_DELTA(vertex_element.GetNode(0)->GetPoint()[0], 0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(0)->GetPoint()[1], 0.5*sqrt(3.0), 1e-9);

        TS_ASSERT_DELTA(vertex_element.GetNode(1)->GetPoint()[0], -0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(1)->GetPoint()[1], 0.5*sqrt(3.0), 1e-9);

        TS_ASSERT_DELTA(vertex_element.GetNode(2)->GetPoint()[0], -0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(2)->GetPoint()[1], -0.5*sqrt(3.0), 1e-9);

        TS_ASSERT_DELTA(vertex_element.GetNode(3)->GetPoint()[0], 0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(3)->GetPoint()[1], -0.5*sqrt(3.0), 1e-9);

        // Add new node
        Node<2>* p_new_node = new Node<2>(4, false, 0.0, 0.0);
        vertex_element.AddNode(3, p_new_node); // Add node at (0,0) between nodes 3 and 0

        // Test node is added
        TS_ASSERT_EQUALS(vertex_element.GetNumNodes(), 5u);

        // Test other nodes are updated
        TS_ASSERT_DELTA(vertex_element.GetNode(0)->GetPoint()[0], 0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(0)->GetPoint()[1], 0.5*sqrt(3.0), 1e-9);

        TS_ASSERT_DELTA(vertex_element.GetNode(1)->GetPoint()[0], -0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(1)->GetPoint()[1], 0.5*sqrt(3.0), 1e-9);

        TS_ASSERT_DELTA(vertex_element.GetNode(2)->GetPoint()[0], -0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(2)->GetPoint()[1], -0.5*sqrt(3.0), 1e-9);

        TS_ASSERT_DELTA(vertex_element.GetNode(3)->GetPoint()[0], 0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(3)->GetPoint()[1], -0.5*sqrt(3.0), 1e-9);

        TS_ASSERT_DELTA(vertex_element.GetNode(4)->GetPoint()[0], 0.0, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(4)->GetPoint()[1], 0.0, 1e-9);

        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
        delete p_new_node;
    }

    void TestMarkAsDeleted()
    {
        // Create nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));

        // Create element
        VertexElement<2,2> vertex_element(0, nodes);
        vertex_element.RegisterWithNodes();

        for (unsigned i=0; i<nodes.size(); i++)
        {
            TS_ASSERT_EQUALS(vertex_element.GetNode(i)->GetNumContainingElements(), 1u);
        }

        vertex_element.MarkAsDeleted();

        for (unsigned i=0; i<nodes.size(); i++)
        {
            TS_ASSERT_EQUALS(vertex_element.GetNode(i)->GetNumContainingElements(), 0u);
        }

        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestUpdateNode()
    {
        // Create nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));

        // Create element
        VertexElement<2,2> vertex_element(0, nodes);
        vertex_element.RegisterWithNodes();

        TS_ASSERT_DELTA(vertex_element.GetNode(2)->rGetLocation()[0], 1.0, 1e-12);
        TS_ASSERT_DELTA(vertex_element.GetNode(2)->rGetLocation()[0], 1.0, 1e-12);

        // Update location of node 2
        Node<2>* p_node = new Node<2>(4, false, 1.2, 1.3);
        vertex_element.UpdateNode(2, p_node);

        TS_ASSERT_DELTA(vertex_element.GetNode(2)->rGetLocation()[0], 1.2, 1e-12);
        TS_ASSERT_DELTA(vertex_element.GetNode(2)->rGetLocation()[1], 1.3, 1e-12);

        // Tidy up
        for (unsigned i=0; i<nodes.size(); ++i)
        {
            delete nodes[i];
        }
        delete p_node;
    }

//     void xTestAnticlockwisenessOfNodes() throw(Exception)
//     {
//        // Tests to check that the nodes are anticlockwise when we create element
//        std::vector<Node<2>*> corner_nodes;
//        corner_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
//        corner_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
//        corner_nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
//        corner_nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
//
//        VertexElement<2,2> vertex_element(INDEX_IS_NOT_USED, corner_nodes);
//
//        std::vector<Node<2>*> corner_nodes2;
//        corner_nodes2.push_back(new Node<2>(0, false, 0.0, 0.0));
//        corner_nodes2.push_back(new Node<2>(2, false, 1.0, 1.0));
//        corner_nodes2.push_back(new Node<2>(1, false, 1.0, 0.0));
//        corner_nodes2.push_back(new Node<2>(3, false, 0.0, 1.0));
//
//        VertexElement<2,2> vertex_element2(INDEX_IS_NOT_USED, corner_nodes2);
//
//        for (unsigned i=0; i<corner_nodes.size(); ++i)
//        {
//            delete corner_nodes[i];
//        }
//        for (unsigned i=0; i<corner_nodes2.size(); ++i)
//        {
//            delete corner_nodes2[i];
//        }
//     }

    void TestGetNodeLocalIndex()
    {
        // Create nodes
        std::vector<Node<2>*> nodes;

        // This is a square
        nodes.push_back(new Node<2>(3, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(0, false, 0.0, 1.0));

        // Create element
        VertexElement<2,2> vertex_element(INDEX_IS_NOT_USED, nodes);

        TS_ASSERT_EQUALS(vertex_element.GetNodeLocalIndex(0), 3u);
        TS_ASSERT_EQUALS(vertex_element.GetNodeLocalIndex(1), 2u);
        TS_ASSERT_EQUALS(vertex_element.GetNodeLocalIndex(2), 1u);
        TS_ASSERT_EQUALS(vertex_element.GetNodeLocalIndex(3), 0u);

        vertex_element.DeleteNode(3); // Removes (1,1) node

        TS_ASSERT_EQUALS(vertex_element.GetNodeLocalIndex(0), UINT_MAX);

        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

};
#endif /*TESTVERTEXELEMENT_HPP_*/
