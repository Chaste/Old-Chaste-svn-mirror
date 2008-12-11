/*

Copyright (C) University of Oxford, 2008

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
#include "VertexElement.hpp"
#include "VertexMesh.hpp"
#include "VertexMeshWriter2d.hpp"

class TestVertexMesh : public CxxTest::TestSuite
{
public:
    void TestBasicVertexMesh() throw(Exception)
    {   // Make four nodes to assign to two elements
        std::vector<Node<2>*> basic_nodes;
        basic_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        basic_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        basic_nodes.push_back(new Node<2>(2, false, 1.5, 1.0));
        basic_nodes.push_back(new Node<2>(3, false, 1.0, 2.0));
        basic_nodes.push_back(new Node<2>(4, false, 0.0, 1.0));
        basic_nodes.push_back(new Node<2>(5, false, 2.0, 0.0));
        basic_nodes.push_back(new Node<2>(6, false, 2.0, 3.0));
        
        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;
        
        // Make two triangular elements out of these nodes
        nodes_elem_0.push_back(basic_nodes[0]);
        nodes_elem_0.push_back(basic_nodes[1]);
        nodes_elem_0.push_back(basic_nodes[2]);
        nodes_elem_0.push_back(basic_nodes[3]);
        nodes_elem_0.push_back(basic_nodes[4]);
        
        nodes_elem_1.push_back(basic_nodes[2]);
        nodes_elem_1.push_back(basic_nodes[5]);
        nodes_elem_1.push_back(basic_nodes[6]);
        
        std::vector<VertexElement<2,2>*> basic_vertex_elements;
        basic_vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        basic_vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));
        // Make a vertex mesh
        VertexMesh<2,2> basic_vertex_mesh(basic_nodes, basic_vertex_elements);
               
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetNumNodes(), 7u);
        
        TS_ASSERT_DELTA(basic_vertex_mesh.GetNode(2)->rGetLocation()[0],1.5,1e-3);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(1)->GetNode(2)->GetIndex(),6u);
        
        // Check that the nodes know which elements they are in
        std::set<unsigned> temp_list1;
        temp_list1.insert(0u);
        
        // Nodes 1 & 4 only in element 0
        TS_ASSERT_EQUALS(basic_nodes[1]->rGetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(basic_nodes[4]->rGetContainingElementIndices(), temp_list1);
        
        // Node 2 in elements 0 and 1
        temp_list1.insert(1u);
        TS_ASSERT_EQUALS(basic_nodes[2]->rGetContainingElementIndices(), temp_list1);
        
        // Node 5 only in element 1
        std::set<unsigned> temp_list2;
        temp_list2.insert(1u);
        TS_ASSERT_EQUALS(basic_nodes[5]->rGetContainingElementIndices(), temp_list2);
        
        for(unsigned i=0; i<basic_nodes.size(); i++)
        {
            delete basic_nodes[i];
        }

        for(unsigned i=0; i<basic_vertex_elements.size(); i++)
        {
            delete basic_vertex_elements[i];
        }
    }
      
    void TestVertexMeshGenerator() throw(Exception)
    {
        VertexMesh<2,2> mesh(3,3);
    
        VertexMeshWriter2d vertex_mesh_writer("TestVertexMeshGeneration","mesh");
        vertex_mesh_writer.WriteFiles(mesh);
    
        /// todo: write test
    
    }

};    

#endif /*TESTVERTEXMESH_HPP_*/
