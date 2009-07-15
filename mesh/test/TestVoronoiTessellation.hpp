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


#ifndef TESTVORONOITESSELLATION_HPP_
#define TESTVORONOITESSELLATION_HPP_

#include "UblasCustomFunctions.hpp"
#include <cxxtest/TestSuite.h>
#include "VoronoiCell.hpp"
#include "VoronoiTessellation.hpp"
#include "MutableMesh.hpp"
#include "Exception.hpp"
#include "TrianglesMeshWriter.hpp"

#include <cmath>
#include <vector>

class TestVoronoiTessellation : public CxxTest::TestSuite
{
public:

    void TestGenerateVerticesFromElementCircumcentres() throw (Exception)
    {
        // Create mutable tetrahedral mesh which is Delauny
        std::vector<Node<3> *> nodes;
        nodes.push_back(new Node<3>(0, true,  1.0,  1.0,  1.0));
        nodes.push_back(new Node<3>(1, true, -1.0, -1.0,  1.0));
        nodes.push_back(new Node<3>(2, true, -1.0,  1.0, -1.0));
        nodes.push_back(new Node<3>(3, true,  1.0, -1.0, -1.0));
        nodes.push_back(new Node<3>(4, false, 0.0,  0.0,  0.0));

        MutableMesh<3,3> mesh(nodes);

        // Note that the Voronois tessellation is not unique for this
        // mesh since 4 points are co-spherical.  We need to check
        // how the mesher is breaking ties.
        Element<3,3> *p_element = mesh.GetElement(0);
        //TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0), 3u);//Older tetgen
        //TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1), 0u);//Older tetgen
        //TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2), 2u);//Older tetgen
        //TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(3), 4u);//Older tetgen
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0), 4u);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1), 1u);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2), 0u);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(3), 2u);

        // Create Voronoi Tessellation
        VoronoiTessellation<3> tessellation(mesh);

        tessellation.GenerateVerticesFromElementCircumcentres();

        TS_ASSERT_EQUALS(tessellation.GetNumVertices(), 8u);

        c_vector<double,3> this_vertex = *(tessellation.GetVertex(2));

        TS_ASSERT_DELTA(this_vertex[0],  1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[1],  1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[2], -1.5, 1e-7);


        this_vertex = *(tessellation.mVertices[2]);

        TS_ASSERT_DELTA(this_vertex[0],  1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[1],  1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[2], -1.5, 1e-7);

        this_vertex = *(tessellation.mVertices[3]);

        TS_ASSERT_DELTA(this_vertex[0],  1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[1], -1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[2],  1.5, 1e-7);

        this_vertex = *(tessellation.mVertices[0]);

        TS_ASSERT_DELTA(this_vertex[0], -1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[1],  1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[2],  1.5, 1e-7);

        this_vertex = *(tessellation.mVertices[1]);

        TS_ASSERT_DELTA(this_vertex[0], -1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[1], -1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[2], -1.5, 1e-7);
    }

    void TestSimpleTessellation() throw (Exception)
    {
        // Create mutable tetrahedral mesh which is Delauny
        std::vector<Node<3> *> nodes;

        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true,  1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(2, true,  1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(3, true,  0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(4, false, 0.5, 0.5, 0.5));

        MutableMesh<3,3> mesh(nodes);
        //TrianglesMeshWriter<3,3> mesh_writer("","Simple_Delaunay_Tet");
        //mesh_writer.WriteFilesUsingMesh(mesh);
        TS_ASSERT(mesh.CheckVoronoi());

        // Create expected cell
        c_vector<double, 3> vertex1;
        vertex1(0) = -0.2500;
        vertex1(1) = -0.2500;
        vertex1(2) = 1.2500;
        c_vector<double, 3> vertex2;
        vertex2(0) = 1.2500;
        vertex2(1) = -0.2500;
        vertex2(2) = -0.2500;
        c_vector<double, 3> vertex3;
        vertex3(0) = -0.2500;
        vertex3(1) = 1.2500;
        vertex3(2) = -0.2500;
        c_vector<double, 3> vertex4;
        vertex4(0) = 1.2500;
        vertex4(1) = 1.2500;
        vertex4(2) = 1.2500;

        Face<3> face1;
        face1.AddVertex(&vertex2);
        face1.AddVertex(&vertex3);
        face1.AddVertex(&vertex4);
        Face<3> face2;
        face2.AddVertex(&vertex1);
        face2.AddVertex(&vertex4);
        face2.AddVertex(&vertex3);
        Face<3> face3;
        face3.AddVertex(&vertex1);
        face3.AddVertex(&vertex2);
        face3.AddVertex(&vertex4);
        Face<3> face4;
        face4.AddVertex(&vertex1);
        face4.AddVertex(&vertex3);
        face4.AddVertex(&vertex2);

        VoronoiCell cell;
        cell.AddFace(&face1);
        cell.AddOrientation(true);
        cell.AddFace(&face2);
        cell.AddOrientation(true);
        cell.AddFace(&face3);
        cell.AddOrientation(true);
        cell.AddFace(&face4);
        cell.AddOrientation(true);

        // Create Voronoi Tesselation
        VoronoiTessellation<3> tessellation(mesh);

        // Check tesellation is correct
        for (unsigned cell_index=0; cell_index<mesh.GetNumNodes(); cell_index++)
        {
            if ( mesh.GetNode(cell_index)->IsBoundaryNode() )
            {
                // For a boundary node we get a null cell
                TS_ASSERT_EQUALS(tessellation.rGetCell(cell_index).GetNumFaces(), 0u);
            }
            else
            {
                // Only node 4 is a non-boundary node
                TS_ASSERT_EQUALS(cell_index, 4u);

                // This gives the expected cell
                TS_ASSERT_EQUALS(tessellation.rGetCell(cell_index), cell);

                VoronoiCell cell = tessellation.rGetCell(cell_index);

                // There will only be one non-boundary cell so we know what the position is.
                TS_ASSERT_DELTA(cell.rGetVoronoiCellCentre()[0], 0.5, 1e-5);
                TS_ASSERT_DELTA(cell.rGetVoronoiCellCentre()[1], 0.5, 1e-5);
                TS_ASSERT_DELTA(cell.rGetVoronoiCellCentre()[2], 0.5, 1e-5);
            }
        }

        // Coverage of GetNumCells()
        TS_ASSERT_EQUALS(tessellation.GetNumCells(), mesh.GetNumNodes());
    }

    void TestTessellation2d() throw (Exception)
    {
        std::vector<Node<2> *> nodes;
        nodes.push_back(new Node<2>(0, true, 0,   0));
        nodes.push_back(new Node<2>(0, true, 0,   1));
        nodes.push_back(new Node<2>(0, true, 1,   1));
        nodes.push_back(new Node<2>(0, true, 1,   0));
        nodes.push_back(new Node<2>(0, true, 0.5, 0.5));

        MutableMesh<2,2> mesh(nodes);

        TS_ASSERT(mesh.CheckVoronoi());

        // Create Voronoi Tesselation
        VoronoiTessellation<2> tessellation(mesh);

        Face<2> expected_face;
        c_vector<double, 2> vertex1;
        vertex1(0) = 0.5;
        vertex1(1) = 0.0;
        expected_face.AddVertex(&vertex1);
        c_vector<double, 2> vertex2;
        vertex2(0) = 1.0;
        vertex2(1) = 0.5;
        expected_face.AddVertex(&vertex2);
        c_vector<double, 2> vertex3;
        vertex3(0) = 0.5;
        vertex3(1) = 1.0;
        expected_face.AddVertex(&vertex3);
        c_vector<double, 2> vertex4;
        vertex4(0) = 0.0;
        vertex4(1) = 0.5;
        expected_face.AddVertex(&vertex4);

        TS_ASSERT_EQUALS(tessellation.rGetFace(4), expected_face);
        TS_ASSERT_EQUALS(tessellation.rGetFace(4).GetNumVertices(), 4u);

        std::vector< c_vector<double, 2>*> vertices_of_face4 = tessellation.rGetFace(4).GetVertices();
        c_vector<double, 2> first_vertex_of_face4 = *(vertices_of_face4[0]);
        TS_ASSERT_DELTA( first_vertex_of_face4(0), 0.5, 1e-4);

        // Calculate length of voronoi edge between nodes 4 and 2
        TS_ASSERT_DELTA(tessellation.GetEdgeLength(4u, 2u), pow(2.0, -0.5), 1e-7);
        TS_ASSERT_EQUALS(tessellation.GetNumFaces(), 5u);
    }

    void TestTessellation2dComplex() throw (Exception)
    {
        std::vector<Node<2> *> nodes;
        nodes.push_back(new Node<2>(0, true,  0.0, 0));
        nodes.push_back(new Node<2>(0, true,  0.0, 1));
        nodes.push_back(new Node<2>(0, true, -1.0, 0));
        nodes.push_back(new Node<2>(0, true,  1.0, 0));
        nodes.push_back(new Node<2>(0, true,  0.5, -pow(3,0.5)/2.0));
        nodes.push_back(new Node<2>(0, true, -0.5, -pow(3,0.5)/2.0));

        MutableMesh<2,2> mesh(nodes);

        TS_ASSERT(mesh.CheckVoronoi());

        // Create Voronoi tessellation
        VoronoiTessellation<2> tessellation(mesh);

        TS_ASSERT_DELTA(tessellation.GetEdgeLength(0u, 1u), 1.0, 1e-6);
        TS_ASSERT_DELTA(tessellation.GetEdgeLength(0u, 2u), 0.5 + pow(3,-0.5)/2.0, 1e-6);
        TS_ASSERT_DELTA(tessellation.GetEdgeLength(0u, 3u), 0.5 + pow(3,-0.5)/2.0, 1e-6);
        TS_ASSERT_DELTA(tessellation.GetEdgeLength(0u, 4u), pow(3,-0.5), 1e-6);
        TS_ASSERT_DELTA(tessellation.GetEdgeLength(0u, 5u), pow(3,-0.5), 1e-6);

        TS_ASSERT_DELTA(tessellation.rGetFace(0).GetArea(), pow(3, 0.5)/4.0+0.5, 1e-6);
        TS_ASSERT_DELTA(tessellation.rGetFace(0).GetPerimeter(), 2.0 + pow(3, 0.5), 1e-6);
        TS_ASSERT_DELTA(tessellation.GetFaceArea(0),  pow(3, 0.5)/4.0+0.5, 1e-6);
        TS_ASSERT_DELTA(tessellation.GetFacePerimeter(0), 2.0 + pow(3, 0.5), 1e-6);
   
    }

    void TestOrderVerticesAntiClockwise()
    {
        std::vector<Node<2> *> nodes;
        nodes.push_back(new Node<2>(0, true,  0,  0));
        nodes.push_back(new Node<2>(0, true,  0,  1));
        nodes.push_back(new Node<2>(0, true,  1,  1));
        nodes.push_back(new Node<2>(0, true,  1,  0));
        nodes.push_back(new Node<2>(0, true, 0.5, 0.5));

        MutableMesh<2,2> mesh(nodes);

        TS_ASSERT(mesh.CheckVoronoi());

        // Create Voronoi tessellation
        VoronoiTessellation<2> tessellation(mesh);

        Face<2> original_face = tessellation.rGetFace(0);
        Face<2> *p_reordered_face = tessellation.mFaces[0];

        c_vector<double,2> *p_temp = &(p_reordered_face->rGetVertex(0));

        //std::cout<< "\n Before " << *(p_reordered_face->mVertices[0]);
        //std::cout<< "\n"<< *(original_face.mVertices[0]);
        //std::cout<< "\n Before " << *(p_reordered_face->mVertices[1]);
        //std::cout<< "\n"<< *(original_face.mVertices[1]);

		/*
         * Note that we can only access mVertices directly as 
         * this test class is a friend of VoronoiTessellation
         */
        p_reordered_face->SetVertex(0, &(p_reordered_face->rGetVertex(1)));
        p_reordered_face->SetVertex(1, p_temp);

        //std::cout<< "\n After " << *(p_reordered_face->mVertices[0]);
        //std::cout<< "\n"<< *(original_face.mVertices[0]);
        //std::cout<< "\n After " << *(p_reordered_face->mVertices[1]);
        //std::cout<< "\n"<< *(original_face.mVertices[1]);

        //tessellation.OrderFaceVerticesAntiClockwise(0u);

        TS_ASSERT_EQUALS(*p_reordered_face, original_face);

        //TS_ASSERT_DIFFERS(original_face,original_face);
    }
};

#endif /*TESTVORONOITESSELLATION_HPP_*/
