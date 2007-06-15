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
        
        nodes.push_back(new Node<3>(0, true,  0.0,  0.0,  0.0));
        nodes.push_back(new Node<3>(1, true,  1.0,  1.0,  0.0));
        nodes.push_back(new Node<3>(2, true,  1.0,  0.0,  1.0));
        nodes.push_back(new Node<3>(3, true,  0.0,  1.0,  1.0));
        nodes.push_back(new Node<3>(4, false, 0.5,0.5,0.5));
        
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
        TS_ASSERT_EQUALS(r_voronoi_cells.size(),1u);
        
        VoronoiCell our_voronoi_cell = r_voronoi_cells[0];
        std::set< Node<3>* >& vertices = our_voronoi_cell.GetVertices();

        
        std::vector< c_vector<double, 3> > expected_vertices;
        c_vector<double, 3> vertex;
        vertex(0)= -0.2500;      
        vertex(1)=-0.2500;
        vertex(2)=1.2500;
        expected_vertices.push_back(vertex); 
        
        vertex(0)=1.2500;
        vertex(1)=-0.2500;
        vertex(2)=-0.2500;
        expected_vertices.push_back(vertex);

        vertex(0)= -0.2500;      
        vertex(1)=1.2500;
        vertex(2)=-0.2500;
        expected_vertices.push_back(vertex);
        
        vertex(0)= 1.2500;      
        vertex(1)=1.2500;
        vertex(2)=1.2500;
        expected_vertices.push_back(vertex);
        
        TS_ASSERT_EQUALS(vertices.size(), expected_vertices.size());
        
        // bit of a messy way of testing that expected_vertices == vertices,
        // ignoring ordering..
        for(unsigned i=0; i<expected_vertices.size(); i++)
        {
            bool found = true;
            for (std::set<Node<3>*>::iterator vertex_iterator = vertices.begin();
                 vertex_iterator!=vertices.end();
                 vertex_iterator++)
            {
                Node<3>* our_node = *vertex_iterator;
    
                // assume the our matches the expected vertex
                found = true;
                
                // check component-wise
                for(unsigned j=0; j<3; j++)
                {
                    if( fabs( our_node->rGetLocation()[j] - expected_vertices[i](j)) > 1e-8)
                    {
                        found = false;
                    }
                }
                
                // if the node does match the vertex, break
                if(found) break;
            }
            TS_ASSERT(found);
        }
            
          // create map from expected vertex index to actual vertex index
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
    
    void TestReturnPolarAngle() throw (Exception)
    {
          // Create conforming tetrahedral mesh which is Delaunay
        std::vector<Node<3> *> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0,0.0,0.0));
                
        ConformingTetrahedralMesh<3,3> mesh(nodes);
        
        // Create Voronoi Tesselation
        VoronoiTessellator tessellator(mesh); 
        
        // Four cases to test:
        // x> 0, y>0
        double angle = tessellator.ReturnPolarAngle(1.0,sqrt(3.0));
        TS_ASSERT_DELTA(angle, M_PI/3.0, 1e-7);
        
        angle = tessellator.ReturnPolarAngle(1.0,-sqrt(3.0));
        TS_ASSERT_DELTA(angle, -M_PI/3.0, 1e-7);
        
        angle = tessellator.ReturnPolarAngle(-1.0,sqrt(3.0));
        TS_ASSERT_DELTA(angle, 2.0*M_PI/3.0, 1e-7);
        
        angle = tessellator.ReturnPolarAngle(-1.0,-sqrt(3.0));
        TS_ASSERT_DELTA(angle, -2.0*M_PI/3.0, 1e-7);
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
        
        // These are deleted in the ReMesh in Constructer so don't need 
        // to be deleted here - or do they?
                
        ConformingTetrahedralMesh<3,3> mesh(nodes);
        
        // Create Voronoi Tesselation
        VoronoiTessellator tessellator(mesh);
        
        tessellator.GenerateVerticesFromElementCircumcentres();
        
        std::vector< Node<3>* > vertices = tessellator.rGetVoronoiVertices();
        
        c_vector<double,3> this_vertex = vertices[0]->rGetLocation();
        
        TS_ASSERT_DELTA(this_vertex[0], 1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[1], 1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[2], -1.5, 1e-7);
        
        this_vertex = vertices[1]->rGetLocation();
        
        TS_ASSERT_DELTA(this_vertex[0], 1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[1], -1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[2], 1.5, 1e-7);
        
        this_vertex = vertices[2]->rGetLocation();
        
        TS_ASSERT_DELTA(this_vertex[0], -1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[1], 1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[2], 1.5, 1e-7);
        
        this_vertex = vertices[3]->rGetLocation();
        
        TS_ASSERT_DELTA(this_vertex[0], -1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[1], -1.5, 1e-7);
        TS_ASSERT_DELTA(this_vertex[2], -1.5, 1e-7);
    }
};


#endif /*TESTVORONOITESSELLATOR_HPP_*/
