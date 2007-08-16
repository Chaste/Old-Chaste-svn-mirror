#ifndef TESTVORONOITESSELLATION_HPP_
#define TESTVORONOITESSELLATION_HPP_


#include "UblasCustomFunctions.hpp"
#include <cxxtest/TestSuite.h>
#include "VoronoiCell.hpp"
#include "VoronoiTessellation.cpp"
#include "ConformingTetrahedralMesh.cpp"
#include "Exception.hpp"

#include <cmath>
#include <vector>

class TestVoronoiTessellation : public CxxTest::TestSuite
{
public:
    void TestReturnPolarAngle() throw (Exception)
    {
        // Create conforming tetrahedral mesh which is Delauny
        std::vector<Node<3> *> nodes;
        nodes.push_back(new Node<3>(0, true,  1.0,  1.0,  1.0));
        nodes.push_back(new Node<3>(1, true, -1.0, -1.0,  1.0));
        nodes.push_back(new Node<3>(2, true, -1.0,  1.0, -1.0));
        nodes.push_back(new Node<3>(3, true,  1.0, -1.0, -1.0));
        nodes.push_back(new Node<3>(4, false, 0.0,0.0,0.0));
                
        ConformingTetrahedralMesh<3,3> mesh(nodes);
        
        // Create Voronoi Tesselation
        VoronoiTessellation<3> tessellation(mesh); 
        
        // Four cases to test:
        // x> 0, y>0
        double angle = tessellation.ReturnPolarAngle(1.0,sqrt(3.0));
        TS_ASSERT_DELTA(angle, M_PI/3.0, 1e-7);
        
        angle = tessellation.ReturnPolarAngle(1.0,-sqrt(3.0));
        TS_ASSERT_DELTA(angle, -M_PI/3.0, 1e-7);
        
        angle = tessellation.ReturnPolarAngle(-1.0,sqrt(3.0));
        TS_ASSERT_DELTA(angle, 2.0*M_PI/3.0, 1e-7);
        
        angle = tessellation.ReturnPolarAngle(-1.0,-sqrt(3.0));
        TS_ASSERT_DELTA(angle, -2.0*M_PI/3.0, 1e-7);
        
        // check boundary cases
        angle = tessellation.ReturnPolarAngle(1.0,0.0);
        TS_ASSERT_DELTA(angle, 0.0, 1e-7);
        
        angle = tessellation.ReturnPolarAngle(0.0,1.0);
        TS_ASSERT_DELTA(angle, M_PI/2.0, 1e-7);
        
        angle = tessellation.ReturnPolarAngle(-1.0,0.0);
        TS_ASSERT_DELTA(angle, M_PI, 1e-7);
        
        angle = tessellation.ReturnPolarAngle(0.0,-1.0);
        TS_ASSERT_DELTA(angle, -M_PI/2.0, 1e-7);
        
    }
    
    void TestGenerateVerticesFromElementCircumcentres() throw (Exception)
    {
        // Create conforming tetrahedral mesh which is Delauny
        std::vector<Node<3> *> nodes;
        nodes.push_back(new Node<3>(0, true,  1.0,  1.0,  1.0));
        nodes.push_back(new Node<3>(1, true, -1.0, -1.0,  1.0));
        nodes.push_back(new Node<3>(2, true, -1.0,  1.0, -1.0));
        nodes.push_back(new Node<3>(3, true,  1.0, -1.0, -1.0));
        nodes.push_back(new Node<3>(4, false, 0.0,0.0,0.0));
                
        ConformingTetrahedralMesh<3,3> mesh(nodes);
        
        // Create Voronoi Tesselation
        VoronoiTessellation<3> tessellation(mesh);
        
        tessellation.GenerateVerticesFromElementCircumcentres();
        
        c_vector<double,3> this_vertex = *(tessellation.mVertices[0]);
        
        TS_ASSERT_DELTA(this_vertex[0], 1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[1], 1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[2], -1.5, 1e-7);
        
        this_vertex = *(tessellation.mVertices[1]);
        
        TS_ASSERT_DELTA(this_vertex[0], 1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[1], -1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[2], 1.5, 1e-7);
        
        this_vertex = *(tessellation.mVertices[2]);
        
        TS_ASSERT_DELTA(this_vertex[0], -1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[1], 1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[2], 1.5, 1e-7);
        
        this_vertex = *(tessellation.mVertices[3]);
        
        TS_ASSERT_DELTA(this_vertex[0], -1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[1], -1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[2], -1.5, 1e-7);
    }
    
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
        //TrianglesMeshWriter<3,3> mesh_writer("","Simple_Delaunay_Tet");
        //mesh_writer.WriteFilesUsingMesh(mesh);
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
        
        Face<3> face1;
        face1.mVertices.push_back(&vertex2);
        face1.mVertices.push_back(&vertex3);
        face1.mVertices.push_back(&vertex4);
        Face<3> face2;
        face2.mVertices.push_back(&vertex1);
        face2.mVertices.push_back(&vertex4);
        face2.mVertices.push_back(&vertex3);
        Face<3> face3;
        face3.mVertices.push_back(&vertex1);
        face3.mVertices.push_back(&vertex2);
        face3.mVertices.push_back(&vertex4);
        Face<3> face4;
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
        VoronoiTessellation<3> tessellation(mesh);
        
        // check tesellation is correct
        for (unsigned cell_index=0; cell_index<mesh.GetNumNodes(); cell_index++)
        {
            if ( mesh.GetNode(cell_index)->IsBoundaryNode() )
            {
                // for a boundary node we get a null cell
                TS_ASSERT(tessellation.rGetCell(cell_index).mFaces.size()==0);
            }
            else
            {
                // only node 4 is a non-boundary node
                TS_ASSERT_EQUALS(cell_index, 4u);
                // this gives the expected cell
                TS_ASSERT_EQUALS(tessellation.rGetCell(cell_index), cell);
                
                VoronoiCell cell = tessellation.rGetCell(cell_index);
                // there will only be one non-boundary cell so we know what the position is.
                TS_ASSERT_DELTA(cell.rGetVoronoiCellCentre()[0], 0.5,1e-5);
                TS_ASSERT_DELTA(cell.rGetVoronoiCellCentre()[1], 0.5,1e-5);
                TS_ASSERT_DELTA(cell.rGetVoronoiCellCentre()[2], 0.5,1e-5);
            }
        }
    }
    
    
    void TestTessellation2d() throw (Exception)
    {
        std::vector<Node<2> *> nodes;
        nodes.push_back(new Node<2>(0, true,  0,  0));
        nodes.push_back(new Node<2>(0, true,  0,  1));
        nodes.push_back(new Node<2>(0, true,  1,  1));
        nodes.push_back(new Node<2>(0, true,  1,  0));
        nodes.push_back(new Node<2>(0, true,  0.5,0.5));
                
        ConformingTetrahedralMesh<2,2> mesh(nodes);
                
        TS_ASSERT(mesh.CheckVoronoi());
        
        // Create Voronoi Tesselation
        VoronoiTessellation<2> tessellation(mesh);
        
        Face<2> expected_face;
        c_vector<double, 2> vertex1;
        vertex1(0)= 0.5;      
        vertex1(1)= 0.0;
        expected_face.mVertices.push_back(&vertex1);
        c_vector<double, 2> vertex2;
        vertex2(0)= 1.0;      
        vertex2(1)= 0.5;
        expected_face.mVertices.push_back(&vertex2);
        c_vector<double, 2> vertex3;
        vertex3(0)= 0.5;      
        vertex3(1)= 1.0;
        expected_face.mVertices.push_back(&vertex3);
        c_vector<double, 2> vertex4;
        vertex4(0)= 0.0;      
        vertex4(1)= 0.5;
        expected_face.mVertices.push_back(&vertex4);
        
        TS_ASSERT_EQUALS(*(tessellation.GetFace(4)), expected_face);
        
    }
};

#endif /*TESTVORONOITESSELLATION_HPP_*/
