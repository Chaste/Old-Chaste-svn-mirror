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

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "VertexElement.hpp"
#include "VertexMesh.hpp"
#include "VertexMeshWriter.hpp"

class TestVertexMesh : public CxxTest::TestSuite
{
public:

    void TestBasicVertexMesh() throw(Exception)
    {   
        // Make four nodes to assign to two elements
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
    }
      
    void TestVertexMeshGenerator() throw(Exception)
    {
        VertexMesh<2,2> mesh(5,3); // columns then rows
            
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 15u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 46u);
        
        // Test some random nodes are in the correct place
        TS_ASSERT_DELTA(mesh.GetNode(5)->rGetLocation()[0],4.3301,1e-3);
        TS_ASSERT_DELTA(mesh.GetNode(5)->rGetLocation()[1],0.0,1e-3);
               
        TS_ASSERT_DELTA(mesh.GetNode(43)->rGetLocation()[0],1.732,1e-3);
        TS_ASSERT_DELTA(mesh.GetNode(43)->rGetLocation()[1],3.5,1e-3);
        
        // Test random element has correct nodes
        TS_ASSERT_EQUALS(mesh.GetElement(6)->GetNode(0)->GetIndex(),19u);
        TS_ASSERT_EQUALS(mesh.GetElement(6)->GetNode(3)->GetIndex(),32u);
        
        // Check that the nodes know which elements they are in
        std::set<unsigned> temp_list1;
        temp_list1.insert(0u);
        
        // Nodes 0 & 1 only in element 0
        TS_ASSERT_EQUALS(mesh.GetNode(0)->rGetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(mesh.GetNode(1)->rGetContainingElementIndices(), temp_list1);
        
        // Node 13 in elements 0 and 1 and 5
        temp_list1.insert(1u);
        temp_list1.insert(5u);
        
        TS_ASSERT_EQUALS(mesh.GetNode(13)->rGetContainingElementIndices(), temp_list1);
    }
    
    void noTestT1Swap() throw(Exception)
    {
        // Make 6 nodes to assign to two elements
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
        nodes.push_back(new Node<2>(4, false, 0.5, 0.4));
        nodes.push_back(new Node<2>(5, false, 0.5, 0.6));
        
        
        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1, nodes_elem_2, nodes_elem_3;
        
        // Make two triangular and two romboid elements out of these nodes
        nodes_elem_0.push_back(nodes[2]);
        nodes_elem_0.push_back(nodes[3]);
        nodes_elem_0.push_back(nodes[5]);
        
        nodes_elem_1.push_back(nodes[2]);
        nodes_elem_1.push_back(nodes[5]);
        nodes_elem_1.push_back(nodes[4]);
        nodes_elem_1.push_back(nodes[1]);
        
        nodes_elem_0.push_back(nodes[1]);
        nodes_elem_0.push_back(nodes[4]);
        nodes_elem_0.push_back(nodes[0]);
        
        nodes_elem_1.push_back(nodes[0]);
        nodes_elem_1.push_back(nodes[4]);
        nodes_elem_1.push_back(nodes[5]);
        nodes_elem_1.push_back(nodes[3]);
        
        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));
        vertex_elements.push_back(new VertexElement<2,2>(2, nodes_elem_2));
        vertex_elements.push_back(new VertexElement<2,2>(3, nodes_elem_3));

        // Make a vertex mesh
        VertexMesh<2,2> vertex_mesh(nodes, vertex_elements);
               
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u); 
        
        // Perform a T1 swap on nodes 4 and 5
        vertex_mesh.PerformT1Swap(vertex_mesh.GetNode(4), vertex_mesh.GetNode(4));
        
        // Test some random nodes are in the correct place
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0],0.4,1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1],0.5,1e-8);
               
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0],0.6,1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1],0.5,1e-3);
        
        // Test random elements have correct nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(0)->GetIndex(),2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(1)->GetIndex(),3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(2)->GetIndex(),5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(3)->GetIndex(),4u);
        
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(0)->GetIndex(),2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(1)->GetIndex(),4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(2)->GetIndex(),1u);        
    }
    
    /*
     * This tests that a 'dummy' archive function does not throw any errors.
     */
    void TestArchiveVertexMesh()
    {        
        std::string dirname = "archive";
        OutputFileHandler handler(dirname, false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "vertex_mesh_base.arch";
        
        std::string mesh_filename = "vertex_mesh";
        std::string mesh_pathname = handler.GetOutputDirectoryFullPath() + mesh_filename;
      
        VertexMesh<2,2>* const p_mesh = new VertexMesh<2,2>(5,3);
        /*
         * You need the const above to stop a BOOST_STATIC_ASSERTION failure.
         * This is because the serialization library only allows you to save tracked
         * objects while the compiler considers them const, to prevent the objects 
         * changing during the save, and so object tracking leading to wrong results.
         * 
         * E.g. A is saved once via pointer, then changed, then saved again. The second
         * save notes that A was saved before, so doesn't write its data again, and the
         * change is lost.
         */

        // Create an ouput archive
        {
            TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 46u);
            TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 15u);
 
            // Save the mesh data using mesh writers
            VertexMeshWriter<2,2> mesh_writer(dirname, mesh_filename, false);
            mesh_writer.WriteFilesUsingMesh(*p_mesh);

            // Archive the mesh
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // We have to serialize via a pointer here, or the derived class information is lost
            output_arch << p_mesh;
        }

        {
            // De-serialize and compare
            VertexMesh<2,2>* p_mesh2;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_mesh2;
            
            // Re-initialise the mesh
            p_mesh2->Clear();
            VertexMeshReader2d mesh_reader(mesh_pathname);
            p_mesh2->ConstructFromMeshReader(mesh_reader);

            // Compare the loaded mesh against the original

            TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), p_mesh2->GetNumNodes());

            for (unsigned node_index=0; node_index<p_mesh->GetNumNodes(); node_index++)
            {
                Node<2> *p_node = p_mesh->GetNode(node_index);
                Node<2> *p_node2 = p_mesh2->GetNode(node_index);

                TS_ASSERT_EQUALS(p_node->IsDeleted(), p_node2->IsDeleted());
                TS_ASSERT_EQUALS(p_node->GetIndex(), p_node2->GetIndex());
                TS_ASSERT_EQUALS(p_node->IsBoundaryNode(), p_node2->IsBoundaryNode());

                for (unsigned dimension=0; dimension<2; dimension++)
                {
                    TS_ASSERT_DELTA(p_node->rGetLocation()[dimension], p_node2->rGetLocation()[dimension], 1e-4);
                }
            }

            TS_ASSERT_EQUALS(p_mesh->GetNumElements(), p_mesh2->GetNumElements());

            for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
            {
                TS_ASSERT_EQUALS(p_mesh->GetElement(elem_index)->GetNumNodes(),
                                 p_mesh2->GetElement(elem_index)->GetNumNodes());

                for (unsigned local_index=0; local_index<p_mesh->GetElement(elem_index)->GetNumNodes(); local_index++)
                {
                    TS_ASSERT_EQUALS(p_mesh->GetElement(elem_index)->GetNodeGlobalIndex(local_index),
                                     p_mesh2->GetElement(elem_index)->GetNodeGlobalIndex(local_index));
                }
            }

            // Tidy up
            delete p_mesh2;
        }
    }
    
    void TestMeshConstructionFromMeshReader(void)
    {
        VertexMeshReader2d mesh_reader("notforrelease-cancer/test/data/TestVertexMesh/vertex_mesh");
        VertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);

        // Check some node co-ordinates
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(2)->GetPoint()[0], 1.5, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(2)->GetPoint()[1], 1.0, 1e-6);

        // Check first element has the right nodes
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(0), 2u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(1), 5u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(2), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(1), mesh.GetNode(5));
    }

    void TestMeshConstructionFromMeshReaderIndexedFromOne(void)
    {
        VertexMeshReader2d mesh_reader("notforrelease-cancer/test/data/TestVertexMesh/vertex_mesh_elements_indexed_from_1");
        VertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);

        // Check some node co-ordinates
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(2)->GetPoint()[0], 1.5, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(2)->GetPoint()[1], 1.0, 1e-6);

        // Check first element has the right nodes
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(0), 2u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(1), 5u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(2), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(1), mesh.GetNode(5));
    }


    void TestSetNode()
    {
        // Create mesh
        VertexMeshReader2d mesh_reader("notforrelease-cancer/test/data/TestVertexMesh/vertex_mesh");
        VertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        ChastePoint<2> point = mesh.GetNode(3)->GetPoint();
        TS_ASSERT_DELTA(point[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(point[1], 2.0, 1e-6);

        // Nudge node
        point.SetCoordinate(0, 1.1);
        mesh.SetNode(3, point);
        
        ChastePoint<2> point2 = mesh.GetNode(3)->GetPoint();
        TS_ASSERT_DELTA(point2[0], 1.1, 1e-6);
        TS_ASSERT_DELTA(point2[1], 2.0, 1e-6);
        
        // Nudge node again
        point.SetCoordinate(1, 1.9);
        mesh.SetNode(3, point);
        
        ChastePoint<2> point3 = mesh.GetNode(3)->GetPoint();
        TS_ASSERT_DELTA(point3[0], 1.1, 1e-6);
        TS_ASSERT_DELTA(point3[1], 1.9, 1e-6);
    }
    
    
    void TestAddNodeAndReMeshMethods(void)
    {
        // Create mesh
        VertexMeshReader2d mesh_reader("notforrelease-cancer/test/data/TestVertexMesh/vertex_mesh");
        VertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);
        
        // Test the AddNode() method
        c_vector<double, 2> point;
        point[0] = 1.5;
        point[1] = 0.0;        
        Node<2>* p_node = new Node<2>(7u, point);
        unsigned new_index = mesh.AddNode(p_node);

        TS_ASSERT_EQUALS(new_index, 7u);
        TS_ASSERT_DELTA(mesh.GetNode(7u)->rGetLocation()[0], 1.5, 1e-7);
        TS_ASSERT_DELTA(mesh.GetNode(7u)->rGetLocation()[1], 0.0, 1e-7);

        // Test the AddNode() and Remesh() methods
        point[0] = 2.5;
        point[1] = 1.5;
        Node<2>* p_node2 = new Node<2>(8u, point);

        NodeMap map(mesh.GetNumElements());
        new_index = mesh.AddNode(p_node2);
        mesh.ReMesh(); // call the version of ReMesh() that doesn't use a map

        TS_ASSERT_EQUALS(new_index, 8u);
        TS_ASSERT_DELTA(mesh.GetNode(8u)->rGetLocation()[0], 2.5 ,1e-7);
        TS_ASSERT_DELTA(mesh.GetNode(8u)->rGetLocation()[1], 1.5 ,1e-7);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 9u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);
        
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNumNodes(), 5u);
    }

    
};    

#endif /*TESTVERTEXMESH_HPP_*/
