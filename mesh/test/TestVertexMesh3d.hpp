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

#ifndef TESTVERTEXMESH_HPP_
#define TESTVERTEXMESH_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "VertexMesh3d.hpp"
#include "VertexMeshWriter.hpp"

class TestVertexMesh3d : public CxxTest::TestSuite
{
public:

    void TestBasic3dVertexMesh()
    {
        // Make 8 nodes to assign to a cube and a pyramid element
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, false, 0.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(4, false, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(5, false, 0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(6, false, 1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(7, false, 1.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(8, false, 0.5, 0.5, 1.5));

        std::vector<Node<3>*> nodes_face_0, nodes_face_1, nodes_face_2, nodes_face_3, nodes_face_4, nodes_face_5,
                              nodes_face_6, nodes_face_7, nodes_face_8, nodes_face_9;

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

        // Make 4 triangular faces out of these nodes
        nodes_face_6.push_back(nodes[6]);
        nodes_face_6.push_back(nodes[7]);
        nodes_face_6.push_back(nodes[8]);

        nodes_face_7.push_back(nodes[6]);
        nodes_face_7.push_back(nodes[8]);
        nodes_face_7.push_back(nodes[3]);

        nodes_face_8.push_back(nodes[3]);
        nodes_face_8.push_back(nodes[8]);
        nodes_face_8.push_back(nodes[5]);

        nodes_face_9.push_back(nodes[5]);
        nodes_face_9.push_back(nodes[8]);
        nodes_face_9.push_back(nodes[7]);

        // Make the faces
        std::vector<VertexElement<2,3>*> faces;

        faces.push_back(new VertexElement<2,3>(0, nodes_face_0));
        faces.push_back(new VertexElement<2,3>(1, nodes_face_1));
        faces.push_back(new VertexElement<2,3>(2, nodes_face_2));
        faces.push_back(new VertexElement<2,3>(3, nodes_face_3));
        faces.push_back(new VertexElement<2,3>(4, nodes_face_4));
        faces.push_back(new VertexElement<2,3>(5, nodes_face_5));
        faces.push_back(new VertexElement<2,3>(6, nodes_face_6));
        faces.push_back(new VertexElement<2,3>(7, nodes_face_7));
        faces.push_back(new VertexElement<2,3>(8, nodes_face_8));
        faces.push_back(new VertexElement<2,3>(9, nodes_face_9));

        // Make the elements
        std::vector<VertexElement<2,3>*> faces_element_0, faces_element_1;
        std::vector<bool> orientations_0, orientations_1;

        // Cube element
        faces_element_0.push_back(faces[0]);
        faces_element_0.push_back(faces[1]);
        faces_element_0.push_back(faces[2]);
        faces_element_0.push_back(faces[3]);
        faces_element_0.push_back(faces[4]);
        faces_element_0.push_back(faces[5]);

        orientations_0.push_back(true);
        orientations_0.push_back(true);
        orientations_0.push_back(true);
        orientations_0.push_back(true);
        orientations_0.push_back(true);
        orientations_0.push_back(true);

        // Pyramid element
        faces_element_1.push_back(faces[6]);
        faces_element_1.push_back(faces[7]);
        faces_element_1.push_back(faces[8]);
        faces_element_1.push_back(faces[9]);
        faces_element_1.push_back(faces[5]);

        orientations_1.push_back(true);
        orientations_1.push_back(true);
        orientations_1.push_back(true);
        orientations_1.push_back(true);
        orientations_1.push_back(false);

        std::vector<VertexElement<3,3>*> elements;
        elements.push_back(new VertexElement<3,3>(0, faces_element_0, orientations_0));
        elements.push_back(new VertexElement<3,3>(1, faces_element_1, orientations_1));

        VertexMesh3d mesh(nodes, faces, elements);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 9u);
        TS_ASSERT_EQUALS(mesh.GetNumFaces(), 10u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), 2u);

        // Test location of random node
        TS_ASSERT_DELTA(mesh.GetNode(2)->rGetLocation()[0], 0.0, 1e-3);
        TS_ASSERT_DELTA(mesh.GetNode(2)->rGetLocation()[1], 1.0, 1e-3);
        TS_ASSERT_DELTA(mesh.GetNode(2)->rGetLocation()[2], 0.0, 1e-3);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(2)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(2)->GetIndex(), 6u);

        // Check that the nodes know which elements they are in.
        std::set<unsigned> temp_list1;
        temp_list1.insert(0u);

        // Nodes 0, 1, 2 and 4 are only in element 0.
        TS_ASSERT_EQUALS(nodes[0]->rGetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(nodes[1]->rGetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(nodes[2]->rGetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(nodes[4]->rGetContainingElementIndices(), temp_list1);

        // Node 3, 5, 6 and 7 are in elements 0 and 1.
        temp_list1.insert(1u);
        TS_ASSERT_EQUALS(nodes[3]->rGetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(nodes[5]->rGetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(nodes[6]->rGetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(nodes[7]->rGetContainingElementIndices(), temp_list1);

        // Node 8 is only in element 1
        std::set<unsigned> temp_list2;
        temp_list2.insert(1u);
        TS_ASSERT_EQUALS(nodes[8]->rGetContainingElementIndices(), temp_list2);

        // Coverage
        TS_ASSERT_EQUALS(mesh.SolveNodeMapping(0), 0u);
        TS_ASSERT_EQUALS(mesh.SolveElementMapping(0), 0u);
        TS_ASSERT_EQUALS(mesh.SolveBoundaryElementMapping(0), 0u);
    }

};

#endif /*TESTVERTEXMESH_HPP_*/
