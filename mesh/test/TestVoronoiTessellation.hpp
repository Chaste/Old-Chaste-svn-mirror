#ifndef TESTVORONOITESSELLATION_HPP_
#define TESTVORONOITESSELLATION_HPP_


#include "UblasCustomFunctions.hpp"
#include <cxxtest/TestSuite.h>
#include "VoronoiCell.hpp"
#include "VoronoiTessellation.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include "TrianglesMeshWriter.cpp"
#include "Exception.hpp"

#include <cmath>
#include <vector>

class TestVoronoiTessellation : public CxxTest::TestSuite
{
public:
    void TestSimpleTessellation() throw (Exception)
    {
        
        // Create conforming tetrahedral mesh which is Delauny
        std::vector<Node<3> *> nodes;
        
        nodes.push_back(new Node<3>(0, true,  0.0,  0.0,  0.0));
        nodes.push_back(new Node<3>(1, true,  1.0,  1.0,  0.0));
        nodes.push_back(new Node<3>(2, true,  1.0,  0.0,  1.0));
        nodes.push_back(new Node<3>(3, true,  0.0,  1.0,  1.0));
        nodes.push_back(new Node<3>(4, false, 0.5,0.5,0.5));
        
        ConformingTetrahedralMesh<3,3> mesh(nodes);
        TrianglesMeshWriter<3,3> mesh_writer("","Simple_Delaunay_Tet");
        mesh_writer.WriteFilesUsingMesh(mesh);
        TS_ASSERT(mesh.CheckVoronoi());

        // create expected cell
        c_vector<double, 3> vertex1;
        vertex1(0)= -0.2500;      
        vertex1(1)=-0.2500;
        vertex1(2)=1.2500;
        c_vector<double, 3> vertex2;        
        vertex2(0)=1.2500;
        vertex2(1)=-0.2500;
        vertex2(2)=-0.2500;
        c_vector<double, 3> vertex3; 
        vertex3(0)= -0.2500;      
        vertex3(1)=1.2500;
        vertex3(2)=-0.2500;
        c_vector<double, 3> vertex4;
        vertex4(0)= 1.2500;      
        vertex4(1)=1.2500;
        vertex4(2)=1.2500;
        
        Face face1;
        face1.mVertices.push_back(&vertex2);
        face1.mVertices.push_back(&vertex3);
        face1.mVertices.push_back(&vertex4);
        Face face2;
        face2.mVertices.push_back(&vertex1);
        face2.mVertices.push_back(&vertex4);
        face2.mVertices.push_back(&vertex3);
        Face face3;
        face3.mVertices.push_back(&vertex1);
        face3.mVertices.push_back(&vertex2);
        face3.mVertices.push_back(&vertex4);
        Face face4;
        face4.mVertices.push_back(&vertex1);
        face4.mVertices.push_back(&vertex3);
        face4.mVertices.push_back(&vertex2);
        
        VoronoiCell cell;
        cell.mFaces.push_back(&face1);
        cell.mOrientations.push_back(true);
        cell.mFaces.push_back(&face2);
        cell.mOrientations.push_back(true);
        cell.mFaces.push_back(&face3);
        cell.mOrientations.push_back(true);
        cell.mFaces.push_back(&face4);
        cell.mOrientations.push_back(true);
        
        // Create Voronoi Tesselation
        VoronoiTessellation tessellation(mesh);
       
        for (unsigned cell_index=0; cell_index<mesh.GetNumNodes(); cell_index++)
        {
            if ( mesh.GetNode(cell_index)->IsBoundaryNode() )
            {
                TS_ASSERT(tessellation.GetCell(cell_index).mFaces.size()==0);
            }
            else
            {
                TS_ASSERT_EQUALS(cell_index, 4u);
                TS_ASSERT_EQUALS(tessellation.GetCell(cell_index), cell);
            }
        }
    }
};

#endif /*TESTVORONOITESSELLATION_HPP_*/
