#ifndef TESTVORONOITESSELLATOR_HPP_
#define TESTVORONOITESSELLATOR_HPP_


#include "UblasCustomFunctions.hpp"
#include <cxxtest/TestSuite.h>
#include "VoronoiCell.hpp"
#include "VoronoiTessellator.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include "TrianglesMeshWriter.cpp"
#include "Exception.hpp"

#include <cmath>
#include <vector>

class TestVoronoiTessellator : public CxxTest::TestSuite
{
public:
    void TestSimpleTessellation() throw (Exception)
    {
        
        // Create conforming tetrahedral mesh which is Delauny
        std::vector<Node<3> *> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0,0.0,0.0));
        nodes.push_back(new Node<3>(1, true, 1.0,1.0,0.0));
        nodes.push_back(new Node<3>(2, true, 1.0,0.0,1.0));
        nodes.push_back(new Node<3>(3, true, 0.0,1.0,1.0));
        nodes.push_back(new Node<3>(4, false,0.5, 0.5, 0.5));
        // These are deleted in the ReMesh in Constructer so don't need 
        // to be deleted here - or do they?
        
        
        ConformingTetrahedralMesh<3,3> mesh(nodes);
        TrianglesMeshWriter<3,3> mesh_writer("","Simple_Delaunay_Tet");
        mesh_writer.WriteFilesUsingMesh(mesh);
        
        TS_ASSERT(mesh.CheckVoronoi());
        
        // Create Voronoi Tesselation
        VoronoiTessellator tessellator(mesh);
        
        tessellator.Generate();
        
        
        const std::vector<VoronoiCell>& r_voronoi_cells = tessellator.rGetVoronoiCells();
        
        // Check Tesselation
        TS_ASSERT_EQUALS(r_voronoi_cells.size(),0u);
        
//        std::vector< c_vector<double, 3> >& vertices = voronoi_cell(1).GetVertices();
//        
//        std::vector< c_vector<double, 3> > expected_vertices;
//        c_vector<double, 3> vertex;
//        vertex(0)=1.2500;
//        vertex(1)=-0.2500;
//        vertex(2)=0.2500;
//        expected_vertices.push_back(vertex);
//
//        vertex(0)= -0.2500;      
//        vertex(1)=1.2500;
//        vertex(2)=-0.2500;
//        expected_vertices.push_back(vertex);
//        
//        vertex(0)= -0.2500;      
//        vertex(1)=-0.2500;
//        vertex(2)=1.2500;
//        expected_vertices.push_back(vertex);        
//        
//        vertex(0)= 1.2500;      
//        vertex(1)=1.2500;
//        vertex(2)=1.2500;
//        expected_vertices.push_back(vertex);        
//        
//        // create map from expected vertex index to actual vertex index
//        
//        std::map<unsigned, unsigned> index_map;
//        for (unsigned index_e=0; index_e<3; index_e++)
//        {
//            for (unsigned index_a=0; index_a<3; index_a++)
//            {
//                if expected_vertices(index_e)==vertices(index_a)
//                {
//                    index_map[index_e]=index_a;
//                }
//            }
//        }
//        
//        // check map is a permutation
//        TS_ASSERT_EQUALS(index_map.size(), 4);
//        
//        // check faces : TODO
    }

};


#endif /*TESTVORONOITESSELLATOR_HPP_*/
