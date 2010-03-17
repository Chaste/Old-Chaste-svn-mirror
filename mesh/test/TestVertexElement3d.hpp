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


#ifndef TESTVERTEXELEMENT3D_HPP_
#define TESTVERTEXELEMENT3D_HPP_

#include "UblasCustomFunctions.hpp"
#include <cxxtest/TestSuite.h>
#include "VertexElement3d.hpp"

#include <cmath>
#include <vector>

class TestVertexElement3d : public CxxTest::TestSuite
{
public:
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

        std::vector<bool> orientations;
        orientations.push_back(true);
        orientations.push_back(true);
        orientations.push_back(true);
        orientations.push_back(true);
        orientations.push_back(true);
        orientations.push_back(true);

        //Make a cube element out of these faces
        VertexElement3d element(0,nodes, faces, orientations);

        TS_ASSERT_EQUALS(element.GetNumNodes(),8u);
        TS_ASSERT_EQUALS(element.GetNumFaces(),6u);

        TS_ASSERT_EQUALS(element.GetIndex(),0u);

        //Test the position of some random nodes
        TS_ASSERT_DELTA(element.GetFace(0)->GetNode(0)->rGetLocation()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(element.GetFace(0)->GetNode(0)->rGetLocation()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(element.GetFace(0)->GetNode(0)->rGetLocation()[2], 0.0, 1e-6);

        TS_ASSERT_DELTA(element.GetFace(5)->GetNode(2)->rGetLocation()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(element.GetFace(5)->GetNode(2)->rGetLocation()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(element.GetFace(5)->GetNode(2)->rGetLocation()[2], 1.0, 1e-6);
    }

};

#endif /*TESTVERTEXELEMENT3D_HPP_*/
