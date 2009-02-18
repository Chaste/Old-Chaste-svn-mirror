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

#ifndef TESTVERTEXMESH_HPP_
#define TESTVERTEXMESH_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "VertexMeshWriter.hpp"
#include "VertexMesh.hpp"


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
        VertexMesh<2,2> basic_vertex_mesh(basic_nodes, basic_vertex_elements, 0.15, 1.76);

        TS_ASSERT_EQUALS(basic_vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetNumNodes(), 7u);
        
        TS_ASSERT_DELTA(basic_vertex_mesh.GetNode(2)->rGetLocation()[0], 1.5, 1e-3);
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
        
        // Test Set and Get methods
        TS_ASSERT_DELTA(basic_vertex_mesh.GetCellRearrangementThreshold(), 0.15, 1e-4);
        TS_ASSERT_DELTA(basic_vertex_mesh.GetEdgeDivisionThreshold(), 1.76, 1e-4);
        
        basic_vertex_mesh.SetCellRearrangementThreshold(0.03);
        basic_vertex_mesh.SetEdgeDivisionThreshold(3.0);

        TS_ASSERT_DELTA(basic_vertex_mesh.GetCellRearrangementThreshold(), 0.03, 1e-4);
        TS_ASSERT_DELTA(basic_vertex_mesh.GetEdgeDivisionThreshold(), 3.0, 1e-4);
    }


    void TestGetCentroidOfElement() throw(Exception)
    {
        // Create nodes
        std::vector<Node<2>*> nodes;
        unsigned num_nodes = 6;
        for (unsigned i=0; i<num_nodes; i++)
        {
            double theta = 2.0*M_PI*(double)(i)/(double)(num_nodes); 
            nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
        }

        // Create element
        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        // Create mesh
        VertexMesh<2,2> mesh(nodes, elements, 0.05, 2.0);

        // Test GetCentroidOfElement() method
        c_vector<double, 2> centroid = mesh.GetCentroidOfElement(0);

        TS_ASSERT_DELTA(centroid(0), 0.0, 1e-6);
        TS_ASSERT_DELTA(centroid(1), 0.0, 1e-6);
    }


    void TestGetAreaGradientOfElementAtNode()
    {
        // Create nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));

        // Create element
        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        // Create mesh
        VertexMesh<2,2> mesh(nodes, elements, 0.05, 2.0);
        
        // Test GetAreaGradientOfElementAtNode() method at each node

        VertexElement<2,2>* p_element = mesh.GetElement(0);

        c_vector<double, 2> element_area_gradient = mesh.GetAreaGradientOfElementAtNode(p_element, 0);        
        TS_ASSERT_DELTA(element_area_gradient[0], -0.5, 1e-6);
        TS_ASSERT_DELTA(element_area_gradient[1], -0.5, 1e-6);
        
        element_area_gradient = mesh.GetAreaGradientOfElementAtNode(p_element, 1);        
        TS_ASSERT_DELTA(element_area_gradient[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(element_area_gradient[1], -0.5, 1e-6);
        
        element_area_gradient = mesh.GetAreaGradientOfElementAtNode(p_element, 2);        
        TS_ASSERT_DELTA(element_area_gradient[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(element_area_gradient[1], 0.5, 1e-6);
        
        element_area_gradient = mesh.GetAreaGradientOfElementAtNode(p_element, 3);        
        TS_ASSERT_DELTA(element_area_gradient[0], -0.5, 1e-6);
        TS_ASSERT_DELTA(element_area_gradient[1], 0.5, 1e-6);
    }



    void TestVertexElementAreaAndPerimeter()
    {
        // Create nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));

        // Create element
        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        // Create mesh
        VertexMesh<2,2> mesh(nodes, elements, 0.05, 2.0);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 4u);

        // Check nodes have correct indices
        for (unsigned i=0; i<4; i++)
        {
            TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(i), i);
        }
        
        // Test area and perimeter calculations
        TS_ASSERT_DELTA(mesh.GetAreaOfElement(0), 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetPerimeterOfElement(0), 4.0, 1e-6);
    }


    void TestGetPerimeterGradientAtNode()
    {
        // Create nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));

        // Create element
        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        // Create mesh
        VertexMesh<2,2> mesh(nodes, elements, 0.05, 2.0);

        // Test gradient of area evaluated at each node

        VertexElement<2,2>* p_element = mesh.GetElement(0);

        c_vector<double, 2> element_perimeter_gradient = mesh.GetPerimeterGradientOfElementAtNode(p_element, 0);        
        TS_ASSERT_DELTA(element_perimeter_gradient[0], -1.0, 1e-6);
        TS_ASSERT_DELTA(element_perimeter_gradient[1], -1.0, 1e-6);

        element_perimeter_gradient = mesh.GetPerimeterGradientOfElementAtNode(p_element, 1);        
        TS_ASSERT_DELTA(element_perimeter_gradient[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(element_perimeter_gradient[1], -1.0, 1e-6);

        element_perimeter_gradient = mesh.GetPerimeterGradientOfElementAtNode(p_element, 2);        
        TS_ASSERT_DELTA(element_perimeter_gradient[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(element_perimeter_gradient[1], 1.0, 1e-6);

        element_perimeter_gradient = mesh.GetPerimeterGradientOfElementAtNode(p_element, 3);        
        TS_ASSERT_DELTA(element_perimeter_gradient[0], -1.0, 1e-6);
        TS_ASSERT_DELTA(element_perimeter_gradient[1], 1.0, 1e-6);
    }


    void TestMeshGetWidthAndWidthExtremes(void)
    {
        // Create mesh
        VertexMesh<2,2> mesh(3, 3, 0.01, 2.0);

        // Test GetWidthExtremes() method
        c_vector<double,2> width_extremes = mesh.GetWidthExtremes(0u);
        c_vector<double,2> height_extremes = mesh.GetWidthExtremes(1u);

        TS_ASSERT_DELTA(width_extremes[0], 0.0000, 1e-4);
        TS_ASSERT_DELTA(width_extremes[1], 2.8867, 1e-4);

        TS_ASSERT_DELTA(height_extremes[0], 0.0000, 1e-4);
        TS_ASSERT_DELTA(height_extremes[1], 3.5000, 1e-4);

        // Test GetWidth() method
        double width = mesh.GetWidth(0);
        double height = mesh.GetWidth(1);

        TS_ASSERT_DELTA(width, 2.8867, 1e-4);
        TS_ASSERT_DELTA(height, 3.5000, 1e-4);
    }


    void TestVertexElementAreaAndPerimeterOnCircle()
    {
        // Create nodes
        std::vector<Node<2>*> nodes;
        unsigned num_nodes = 1000;
        for (unsigned i=0; i<num_nodes; i++)
        {
            double theta = 2.0*M_PI*(double)(i)/(double)(num_nodes); 
            nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));   
        }

        // Create element
        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        // Create mesh
        VertexMesh<2,2> mesh(nodes, elements, 0.05, 2.0);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), num_nodes);
        
        //  Check nodes have correct indices
        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(i), i);
        }

        // Test area and perimeter calculations
        TS_ASSERT_DELTA(mesh.GetAreaOfElement(0), M_PI, 1e-4);
        TS_ASSERT_DELTA(mesh.GetPerimeterOfElement(0), 2.0*M_PI, 1e-4);
    }


    void TestVertexMeshGenerator() throw(Exception)
    {
        // Create mesh
        VertexMesh<2,2> mesh(5, 3, 0.01, 2.0);
            
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 15u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 46u);
        
        // Test some random nodes are in the correct place
        TS_ASSERT_DELTA(mesh.GetNode(5)->rGetLocation()[0], 4.3301, 1e-3);
        TS_ASSERT_DELTA(mesh.GetNode(5)->rGetLocation()[1], 0.0, 1e-3);
               
        TS_ASSERT_DELTA(mesh.GetNode(43)->rGetLocation()[0], 1.732, 1e-3);
        TS_ASSERT_DELTA(mesh.GetNode(43)->rGetLocation()[1], 3.5, 1e-3);
        
        // Test random element has correct nodes
        TS_ASSERT_EQUALS(mesh.GetElement(6)->GetNode(0)->GetIndex(), 19u);
        TS_ASSERT_EQUALS(mesh.GetElement(6)->GetNode(3)->GetIndex(), 32u);
        
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
    

    void TestMeshConstructionFromMeshReader(void)
    {
        // Create mesh
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
        // Create mesh
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


    void TestAddNodeAndReMesh() throw (Exception)
    {
        // Create mesh
        VertexMeshReader2d mesh_reader("notforrelease-cancer/test/data/TestVertexMesh/vertex_mesh");
        VertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);

        // Choose a node
        ChastePoint<2> point = mesh.GetNode(3)->GetPoint();
        TS_ASSERT_DELTA(point[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(point[1], 2.0, 1e-6);

        // Create a new node close to this node
        point.SetCoordinate(0, 1.1);
        point.SetCoordinate(1, 2.1);
        Node<2>* p_node = new Node<2>(mesh.GetNumNodes(), point);

        unsigned old_num_nodes = mesh.GetNumNodes();

        // Add this new node to the mesh
        unsigned new_index = mesh.AddNode(p_node);
        TS_ASSERT_EQUALS(new_index, old_num_nodes);

        // Remesh to update correspondences
        mesh.SetEdgeDivisionThreshold(1000); // set high threshold to avoid more nodes appearing in the remesh
        NodeMap map(mesh.GetNumElements());
        mesh.ReMesh(map);

        TS_ASSERT_EQUALS(map.Size(), mesh.GetNumElements());
        TS_ASSERT_EQUALS(map.IsIdentityMap(), true);

        // Check that the mesh is updated
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);

        TS_ASSERT_DELTA(mesh.GetNode(new_index)->rGetLocation()[0], 1.1, 1e-7);
        TS_ASSERT_DELTA(mesh.GetNode(new_index)->rGetLocation()[1], 2.1, 1e-7);

        // Now test AddNode() when mDeletedNodeIndices is populated

        // Label node 3 as deleted
        mesh.mDeletedNodeIndices.push_back(3);

        // Create a new node close to this node
        ChastePoint<2> point2;
        point2.SetCoordinate(0, 0.9);
        point2.SetCoordinate(1, 1.9);
        Node<2>* p_node2 = new Node<2>(mesh.GetNumNodes(), point);

        // Add this new node to the mesh
        new_index = mesh.AddNode(p_node2);
        TS_ASSERT_EQUALS(new_index, 3u);
    }


    void TestAddElement() throw (Exception)
    {
        // Make four nodes to assign to two elements
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.5, 1.0));
        nodes.push_back(new Node<2>(3, false, 1.0, 2.0));
        nodes.push_back(new Node<2>(4, false, 0.0, 1.0));
        nodes.push_back(new Node<2>(5, false, 2.0, 0.0));
        nodes.push_back(new Node<2>(6, false, 2.0, 3.0));
        
        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;
        
        // Make two triangular elements out of these nodes
        nodes_elem_0.push_back(nodes[0]);
        nodes_elem_0.push_back(nodes[1]);
        nodes_elem_0.push_back(nodes[2]);
        nodes_elem_0.push_back(nodes[3]);
        nodes_elem_0.push_back(nodes[4]);
        
        nodes_elem_1.push_back(nodes[2]);
        nodes_elem_1.push_back(nodes[5]);
        nodes_elem_1.push_back(nodes[6]);
        
        std::vector<VertexElement<2,2>*> elements;
        VertexElement<2,2> *p_replaced_vertex_element=new VertexElement<2,2>(0, nodes_elem_0);
        elements.push_back(p_replaced_vertex_element);
        elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));
        
        // Make a vertex mesh
        VertexMesh<2,2> mesh(nodes, elements, 0.05, 2.0);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);

        std::vector<Node<2>*> nodes_elem_2, nodes_elem_3;
        
        // Make two triangular elements out of these nodes
        nodes_elem_2.push_back(nodes[6]);
        nodes_elem_2.push_back(nodes[1]);
        nodes_elem_2.push_back(nodes[2]);

        nodes_elem_3.push_back(nodes[0]);
        nodes_elem_3.push_back(nodes[1]);
        nodes_elem_3.push_back(nodes[2]);
        nodes_elem_3.push_back(nodes[3]);

        // Add a new element to the mesh 
        mesh.AddElement(new VertexElement<2,2>(2, nodes_elem_2));

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 3u);


        // Replace element 0 in the mesh 
        mesh.AddElement(new VertexElement<2,2>(0, nodes_elem_3));

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 3u);
        
        //Tidy up
        delete p_replaced_vertex_element;
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
      
        VertexMesh<2,2>* const p_mesh = new VertexMesh<2,2>(5, 3, 0.01, 2.0);
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
            
            
            TS_ASSERT_DELTA(p_mesh2->GetCellRearrangementThreshold(), 0.01, 1e-4);
            TS_ASSERT_DELTA(p_mesh2->GetEdgeDivisionThreshold(), 2.0, 1e-4);

            // Tidy up
            delete p_mesh;
            delete p_mesh2;
        }
    }

    void TestNodesMerging() throw(Exception)
    {
        // Make 4 nodes to assign to 2 elements
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
        nodes.push_back(new Node<2>(4, false, 0.4, 0.0));
        nodes.push_back(new Node<2>(5, false, 0.6, 0.0));
        nodes.push_back(new Node<2>(6, false, 0.4, 0.4));
        nodes.push_back(new Node<2>(7, false, 0.6, 0.6));
        

        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;
                
        // Make two triangular  elements out of these nodes
        nodes_elem_0.push_back(nodes[0]);
        nodes_elem_0.push_back(nodes[4]);
        nodes_elem_0.push_back(nodes[5]);
        nodes_elem_0.push_back(nodes[1]);
        nodes_elem_0.push_back(nodes[2]);
        nodes_elem_0.push_back(nodes[7]);
        nodes_elem_0.push_back(nodes[6]);
        
        nodes_elem_1.push_back(nodes[0]);
        nodes_elem_1.push_back(nodes[6]);
        nodes_elem_1.push_back(nodes[7]);
        nodes_elem_1.push_back(nodes[2]);
        nodes_elem_1.push_back(nodes[3]);
        
        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));
        
        // Make a vertex mesh
        VertexMesh<2,2> vertex_mesh(nodes, vertex_elements);
               
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 8u); 
        
        // Perform a merge on nodes 4 and 5
        std::set<unsigned> containing_element_indices;
        containing_element_indices.insert(0);
        vertex_mesh.PerformNodeMerge(vertex_mesh.GetNode(4), vertex_mesh.GetNode(5), containing_element_indices);

        // Perform a merge on nodes 6 and 7
        containing_element_indices.insert(1);
        vertex_mesh.PerformNodeMerge(vertex_mesh.GetNode(6), vertex_mesh.GetNode(7), containing_element_indices);
                
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 8u);
        
        // Test merged nodes are in the correct place
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.0, 1e-8);
               
        TS_ASSERT_DELTA(vertex_mesh.GetNode(6)->rGetLocation()[0], 0.5, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(6)->rGetLocation()[1], 0.5, 1e-3);
        
        // Test elements have correct nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(1)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(2)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(3)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(4)->GetIndex(), 6u);
        
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(1)->GetIndex(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(2)->GetIndex(), 2u); 
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(3)->GetIndex(), 3u);
        
        // Test Areas and Perimeters of elements 
        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(0), 0.5, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(0), 2+sqrt(2), 1e-6);
        
        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(1), 0.5,1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(1), 2.0+sqrt(2), 1e-6);
    }

    void TestPerformT1Swap() throw(Exception)
    {
        // Make 6 nodes to assign to four elements
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
        nodes.push_back(new Node<2>(4, false, 0.5, 0.4));
        nodes.push_back(new Node<2>(5, false, 0.5, 0.6));

        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1, nodes_elem_2, nodes_elem_3;
        
        // Make two triangular and two rhomboid elements out of these nodes
        nodes_elem_0.push_back(nodes[2]);
        nodes_elem_0.push_back(nodes[3]);
        nodes_elem_0.push_back(nodes[5]);
        
        nodes_elem_1.push_back(nodes[2]);
        nodes_elem_1.push_back(nodes[5]);
        nodes_elem_1.push_back(nodes[4]);
        nodes_elem_1.push_back(nodes[1]);
        
        nodes_elem_2.push_back(nodes[1]);
        nodes_elem_2.push_back(nodes[4]);
        nodes_elem_2.push_back(nodes[0]);
        
        nodes_elem_3.push_back(nodes[0]);
        nodes_elem_3.push_back(nodes[4]);
        nodes_elem_3.push_back(nodes[5]);
        nodes_elem_3.push_back(nodes[3]);
        
        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));
        vertex_elements.push_back(new VertexElement<2,2>(2, nodes_elem_2));
        vertex_elements.push_back(new VertexElement<2,2>(3, nodes_elem_3));

        // Make a vertex mesh
        VertexMesh<2,2> vertex_mesh(nodes, vertex_elements, 0.1); // threshold distance is 0.1 to ease calculations
               
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u); 
        
        // Test areas and perimeters of elements 
        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(0), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(0), 1.0+0.2*sqrt(41.0), 1e-6);
                
        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(1),0.3,1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(1), 1.2+0.2*sqrt(41.0), 1e-6);
        
        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(2),0.2,1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(2), 1.0+0.2*sqrt(41.0), 1e-6);
        
        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(3),0.3,1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(3), 1.2+0.2*sqrt(41.0), 1e-6);
        
        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(3), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(3), 1.2+0.2*sqrt(41.0), 1e-6);

        // Perform a T1 swap on nodes 4 and 5
        std::set<unsigned> containing_element_indices;
        containing_element_indices.insert(0);
        containing_element_indices.insert(1);
        containing_element_indices.insert(2);
        containing_element_indices.insert(3);

        vertex_mesh.PerformT1Swap(vertex_mesh.GetNode(4), vertex_mesh.GetNode(5), containing_element_indices);
        
        // Test moved nodes are in the correct place
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-8);
               
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.4, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.5, 1e-3);
        
        // Test elements have correct nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(0)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(1)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(2)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(3)->GetIndex(), 4u);
        
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(0)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(1)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(2)->GetIndex(), 1u); 
        
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(0)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(1)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(2)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(3)->GetIndex(), 0u);
        
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(1)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(2)->GetIndex(), 3u);       
        
        // Test areas and perimeters of elements 
        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(0), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(0), 1.2+0.2*sqrt(41.0), 1e-6);
        
        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(1), 0.2,1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(1), 1.0+0.2*sqrt(41.0), 1e-6);
        
        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(2), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(2), 1.2+0.2*sqrt(41.0), 1e-6);
        
        TS_ASSERT_DELTA(vertex_mesh.GetAreaOfElement(3), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetPerimeterOfElement(3), 1.0+0.2*sqrt(41.0), 1e-6); 
    }


    void TestReMesh() throw(Exception)
    {
        // Create mesh
        VertexMeshReader2d mesh_reader("notforrelease-cancer/test/data/TestVertexMesh/test_remesh_mesh");
        VertexMesh<2,2> vertex_mesh;
        vertex_mesh.ConstructFromMeshReader(mesh_reader);
        vertex_mesh.SetCellRearrangementThreshold(0.1);
        vertex_mesh.SetEdgeDivisionThreshold(DBL_MAX);
               
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 12u); 
        
        // Call remesh 
        vertex_mesh.ReMesh();
        
        // Test moved nodes are in the correct place

        // Center T1Swap
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-8);
               
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.4, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.5, 1e-8);

        // Bottom left T1Swap
        TS_ASSERT_DELTA(vertex_mesh.GetNode(0)->rGetLocation()[0], 0.1, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(0)->rGetLocation()[1], 0.0, 1e-8);
               
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[0], -0.1, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[1], 0.0, 1e-8);

        // Test elements have correct nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(0)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(1)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(2)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNode(3)->GetIndex(), 4u);
        
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(0)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(1)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(2)->GetIndex(), 1u); 
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNode(3)->GetIndex(), 9u);
        
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(0)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(1)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(2)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(3)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(4)->GetIndex(), 11u);        
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNode(5)->GetIndex(), 7u);         

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(1)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(3)->GetNode(2)->GetIndex(), 3u);
     
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNode(0)->GetIndex(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNode(1)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNode(2)->GetIndex(), 9u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(4)->GetNode(3)->GetIndex(), 1u);
    }


    void TestDivideEdgeIfTooBig() throw(Exception)
    {
        // Create nodes
        Node<2>* p_node0 = new Node<2>(0, false, 1.0, 1.0);
        Node<2>* p_node1 = new Node<2>(1, false, 2.0, 1.0);
        Node<2>* p_node2 = new Node<2>(2, false, 2.0, 3.0);
        Node<2>* p_node3 = new Node<2>(3, false, 1.0, 2.0);

        std::vector<Node<2>*> nodes;        
        nodes.push_back(p_node0);
        nodes.push_back(p_node1);
        nodes.push_back(p_node2);
        nodes.push_back(p_node3);

        // Create element
        VertexElement<2,2>* p_element = new VertexElement<2,2>(0, nodes);
        std::vector<VertexElement<2,2>* > elements;
        elements.push_back(p_element);

        // Create mesh
        VertexMesh<2,2> mesh(nodes, elements);
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 1u);
        
        // Call remesh 
        mesh.ReMesh();

        // Check that the edge between nodes 1 and 2 has divided
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 1u);
        
        TS_ASSERT_DELTA(mesh.GetNode(4)->rGetLocation()[0], 2.0, 1e-8);
        TS_ASSERT_DELTA(mesh.GetNode(4)->rGetLocation()[1], 2.0, 1e-8);
    }


    void TestNeighbouringNodeMethods() throw(Exception)
    {
        // Create mesh
        VertexMesh<2,2> mesh(2, 2, 0.01, 2.0);
            
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 16u);
        
        // Check we have the correct neighbours for node 6
        std::set<unsigned> neighbours = mesh.GetNeighbouringNodeIndices(6);
        
        std::set<unsigned> expected_neighbours;
        expected_neighbours.insert(3);
        expected_neighbours.insert(9);
        expected_neighbours.insert(5);
        
        TS_ASSERT_EQUALS(neighbours, expected_neighbours);
        
        // Check that the only neighbour not also in element 2 is node 3
        std::set<unsigned> neighbours_not_in_elem2 = mesh.GetNeighbouringNodeNotAlsoInElement(6, 2);
        
        TS_ASSERT_EQUALS(neighbours_not_in_elem2.size(), 1u);
        TS_ASSERT_EQUALS(*(neighbours_not_in_elem2.begin()), 3u);
    }
    
    
    void TestDivideEdge()
    {
        // Create nodes
        Node<2>* p_node0 = new Node<2>(0, false, 1.0, 1.0);
        Node<2>* p_node1 = new Node<2>(1, false, 2.0, 1.0);
        Node<2>* p_node2 = new Node<2>(2, false, 2.0, 2.0);
        Node<2>* p_node3 = new Node<2>(3, false, 1.0, 2.0);
        Node<2>* p_node4 = new Node<2>(4, false, 0.5, 1.5);

        std::vector<Node<2>*> nodes_in_element0;        
        nodes_in_element0.push_back(p_node0);
        nodes_in_element0.push_back(p_node1);
        nodes_in_element0.push_back(p_node2);
        nodes_in_element0.push_back(p_node3);

        std::vector<Node<2>*> nodes_in_element1;
        nodes_in_element1.push_back(p_node3);
        nodes_in_element1.push_back(p_node4);
        nodes_in_element1.push_back(p_node0);

        std::vector<Node<2>*> nodes;
        nodes.push_back(p_node0);
        nodes.push_back(p_node1);
        nodes.push_back(p_node2);
        nodes.push_back(p_node3);
        nodes.push_back(p_node4);

        // Create two elements, which share nodes 2 and 3

        VertexElement<2,2>* p_element0 = new VertexElement<2,2>(0, nodes_in_element0);
        VertexElement<2,2>* p_element1 = new VertexElement<2,2>(1, nodes_in_element1);
        
        TS_ASSERT_EQUALS(p_element0->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(p_element1->GetNumNodes(), 3u);
        
        // Create mesh

        std::vector<VertexElement<2,2>* > elements;
        elements.push_back(p_element0);
        elements.push_back(p_element1);
        
        VertexMesh<2,2> mesh(nodes, elements);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNumNodes(), 3u);

        // Divide the edge joining nodes 0 and 1
        mesh.DivideEdge(mesh.GetNode(0), mesh.GetNode(1));
        
        // Test edge is divided

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNumNodes(), 3u);
        TS_ASSERT_DELTA(mesh.GetAreaOfElement(0), 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetPerimeterOfElement(0), 4.0, 1e-6);
        
        // Test other nodes are updated

        TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(0)->GetPoint()[0], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(0)->GetPoint()[1], 1.0, 1e-9);       
        TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(1)->GetPoint()[0], 1.5, 1e-9);
        TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(1)->GetPoint()[1], 1.0, 1e-9);           
        TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(2)->GetPoint()[0], 2.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(2)->GetPoint()[1], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(3)->GetPoint()[0], 2.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(3)->GetPoint()[1], 2.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(4)->GetPoint()[0], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(4)->GetPoint()[1], 2.0, 1e-9);

        TS_ASSERT_DELTA(mesh.GetElement(1)->GetNode(0)->GetPoint()[0], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetElement(1)->GetNode(0)->GetPoint()[1], 2.0, 1e-9);       
        TS_ASSERT_DELTA(mesh.GetElement(1)->GetNode(1)->GetPoint()[0], 0.5, 1e-9);
        TS_ASSERT_DELTA(mesh.GetElement(1)->GetNode(1)->GetPoint()[1], 1.5, 1e-9);           
        TS_ASSERT_DELTA(mesh.GetElement(1)->GetNode(2)->GetPoint()[0], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetElement(1)->GetNode(2)->GetPoint()[1], 1.0, 1e-9);

        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[1], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(1)->GetPoint()[0], 2.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(1)->GetPoint()[1], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(2)->GetPoint()[0], 2.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(2)->GetPoint()[1], 2.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(3)->GetPoint()[0], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(3)->GetPoint()[1], 2.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(4)->GetPoint()[0], 0.5, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(4)->GetPoint()[1], 1.5, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(5)->GetPoint()[0], 1.5, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(5)->GetPoint()[1], 1.0, 1e-9);

        // Divide the edge joining nodes 3 and 0

        mesh.DivideEdge(mesh.GetNode(3), mesh.GetNode(0));
    
        // Test edge is divided

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_DELTA(mesh.GetAreaOfElement(0), 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetPerimeterOfElement(0), 4.0, 1e-6);
                 
        // Test other nodes are updated

        TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(0)->GetPoint()[0], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(0)->GetPoint()[1], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(1)->GetPoint()[0], 1.5, 1e-9);
        TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(1)->GetPoint()[1], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(2)->GetPoint()[0], 2.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(2)->GetPoint()[1], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(3)->GetPoint()[0], 2.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(3)->GetPoint()[1], 2.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(4)->GetPoint()[0], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(4)->GetPoint()[1], 2.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(5)->GetPoint()[0], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(5)->GetPoint()[1], 1.5, 1e-9);

        TS_ASSERT_DELTA(mesh.GetElement(1)->GetNode(0)->GetPoint()[0], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetElement(1)->GetNode(0)->GetPoint()[1], 2.0, 1e-9);       
        TS_ASSERT_DELTA(mesh.GetElement(1)->GetNode(1)->GetPoint()[0], 0.5, 1e-9);
        TS_ASSERT_DELTA(mesh.GetElement(1)->GetNode(1)->GetPoint()[1], 1.5, 1e-9);           
        TS_ASSERT_DELTA(mesh.GetElement(1)->GetNode(2)->GetPoint()[0], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetElement(1)->GetNode(2)->GetPoint()[1], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetElement(1)->GetNode(3)->GetPoint()[0], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetElement(1)->GetNode(3)->GetPoint()[1], 1.5, 1e-9);

        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[1], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(1)->GetPoint()[0], 2.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(1)->GetPoint()[1], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(2)->GetPoint()[0], 2.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(2)->GetPoint()[1], 2.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(3)->GetPoint()[0], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(3)->GetPoint()[1], 2.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(4)->GetPoint()[0], 0.5, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(4)->GetPoint()[1], 1.5, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(5)->GetPoint()[0], 1.5, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(5)->GetPoint()[1], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(6)->GetPoint()[0], 1.0, 1e-9);
        TS_ASSERT_DELTA(mesh.GetNode(6)->GetPoint()[1], 1.5, 1e-9);
    }


    void TestDivideVertexElementGivenNodes() throw(Exception)
    {
        // Make four nodes
        std::vector<Node<2>*> basic_nodes;
        basic_nodes.push_back(new Node<2>(0, false, 2.0, -1.0));
        basic_nodes.push_back(new Node<2>(1, false, 2.0, 1.0));
        basic_nodes.push_back(new Node<2>(2, false, -2.0, 1.0));
        basic_nodes.push_back(new Node<2>(3, false, -2.0, -1.0));

        std::vector<Node<2>*> nodes_elem;

        // Make one rectangular element out of these nodes
        nodes_elem.push_back(basic_nodes[0]);
        nodes_elem.push_back(basic_nodes[1]);
        nodes_elem.push_back(basic_nodes[2]);
        nodes_elem.push_back(basic_nodes[3]);

        std::vector<VertexElement<2,2>*> basic_vertex_elements;
        basic_vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem));

        // Make a vertex mesh
        VertexMesh<2,2> basic_vertex_mesh(basic_nodes, basic_vertex_elements);

        TS_ASSERT_EQUALS(basic_vertex_mesh.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetNumNodes(), 4u);

        // Divide element using two given nodes
        unsigned new_element_index = basic_vertex_mesh.DivideElement(basic_vertex_mesh.GetElement(0), 2, 0);

        TS_ASSERT_EQUALS(new_element_index, 1u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetNumElements(), 2u);

        // Test elements have correct nodes
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(0)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(0)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(0)->GetNode(1)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(0)->GetNode(2)->GetIndex(), 2u);

        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(1)->GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(1)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(1)->GetNode(1)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(1)->GetNode(2)->GetIndex(), 3u);
    }


    void TestDivideVertexElementAlongShortAxis() throw(Exception)
    {
        // Make four nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 2.0, -1.0));
        nodes.push_back(new Node<2>(1, false, 2.0, 1.0));
        nodes.push_back(new Node<2>(2, false, -2.0, 1.0));
        nodes.push_back(new Node<2>(3, false, -2.0, -1.0));
        nodes.push_back(new Node<2>(4, false, 0.0, 2.0));

        // Make a rectangular element out of nodes 0,1,2,3
        std::vector<Node<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[0]);
        nodes_elem_1.push_back(nodes[1]);
        nodes_elem_1.push_back(nodes[2]);
        nodes_elem_1.push_back(nodes[3]);

        // Make a triangular element out of nodes 1,4,2
        std::vector<Node<2>*> nodes_elem_2;
        nodes_elem_2.push_back(nodes[1]);
        nodes_elem_2.push_back(nodes[4]);
        nodes_elem_2.push_back(nodes[2]);

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_1));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_2));

        // Make a vertex mesh
        VertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 5u);
        
        // Divide element 0 along short axis 
        unsigned new_element_index = vertex_mesh.DivideElement(vertex_mesh.GetElement(0));

        TS_ASSERT_EQUALS(new_element_index, vertex_mesh.GetNumElements()-1);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 7u);

        // Now test the position of new nodes
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 1.0, 1e-8);

        TS_ASSERT_DELTA(vertex_mesh.GetNode(6)->rGetLocation()[0], 0.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(6)->rGetLocation()[1], -1.0, 1e-8);

        // Now test the nodes in each element
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNumNodes(), 4u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(0), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(1), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(2), 3u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNodeGlobalIndex(3), 6u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(0), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(1), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(2), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNodeGlobalIndex(3), 5u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(1), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(2), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(2)->GetNodeGlobalIndex(3), 6u);
        
        // Test ownership of the new nodes
        std::set<unsigned> expected_elements_containing_node_5;
        expected_elements_containing_node_5.insert(0);
        expected_elements_containing_node_5.insert(1);
        expected_elements_containing_node_5.insert(2);
        
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(5)->rGetContainingElementIndices(), expected_elements_containing_node_5);
        
        std::set<unsigned> expected_elements_containing_node_6;
        expected_elements_containing_node_6.insert(0);
        expected_elements_containing_node_6.insert(2);
        
        TS_ASSERT_EQUALS(vertex_mesh.GetNode(6)->rGetContainingElementIndices(), expected_elements_containing_node_6);
    }


    void TestDivideVertexElementWithNonRegularElement() throw(Exception)
    {
        // Make six nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(1, false, 2.0, 1.0));
        nodes.push_back(new Node<2>(2, false, 3.0, 2.0));
        nodes.push_back(new Node<2>(3, false, 3.0, 3.0));
        nodes.push_back(new Node<2>(4, false, 1.0, 2.0));

        std::vector<Node<2>*> nodes_elem;

        // Make one rectangular element out of these nodes
        nodes_elem.push_back(nodes[0]);
        nodes_elem.push_back(nodes[1]);
        nodes_elem.push_back(nodes[2]);
        nodes_elem.push_back(nodes[3]);
        nodes_elem.push_back(nodes[4]);

        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes_elem));

        // Make a vertex mesh
        VertexMesh<2,2> mesh(nodes, elements);

        TS_ASSERT_EQUALS(mesh.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 5u);

        // Divide element using two given nodes
        unsigned new_element_index = mesh.DivideElement(mesh.GetElement(0));

        TS_ASSERT_EQUALS(new_element_index, 1u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 7u);

        // Test elements have correct nodes
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(0)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(1)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(2)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(3)->GetIndex(), 6u);

        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(1)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(2)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(3)->GetIndex(), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(4)->GetIndex(), 4u);
        
        // Test locations of new nodes
        TS_ASSERT_DELTA(mesh.GetNode(5)->rGetLocation()[0], 2.4166, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(5)->rGetLocation()[1], 1.4166, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(6)->rGetLocation()[0], 1.5992, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(6)->rGetLocation()[1], 2.2996, 1e-4);
    }


    /// \todo this should be for a non regular polygon so Ixy != 0 (see #825)
    void TestCalculateMomentOfElement() throw(Exception)
    {
        // Create nodes
        std::vector<Node<2>*> nodes;
        unsigned num_nodes = 6;
        for (unsigned i=0; i<num_nodes; i++)
        {
            double theta = 2.0*M_PI*(double)(i)/(double)(num_nodes); 
            nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));   
        }

        // Create element
        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        // Create mesh
        VertexMesh<2,2> mesh(nodes, elements);

        // Test CalculateMomentOfElement() method
        c_vector<double, 3> moments = mesh.CalculateMomentsOfElement(0);

        TS_ASSERT_DELTA(moments(0), 5*sqrt(3)/16, 1e-6);    // Ixx
        TS_ASSERT_DELTA(moments(1), 5*sqrt(3)/16, 1e-6);    // Iyy
        TS_ASSERT_DELTA(moments(2), 0.0, 1e-6);    // Ixy
    }


    void TestGetShortAxisOfElement() throw(Exception)
    {
        // First test
        
        // Create nodes
        std::vector<Node<2>*> nodes1;
       
        // This is a rectangle, centre (0,0), width 2, length 2, parallel to x axis                       
        nodes1.push_back(new Node<2>(0, false,  2.0,  1.0));
        nodes1.push_back(new Node<2>(1, false, -2.0,  1.0));
        nodes1.push_back(new Node<2>(2, false, -2.0, -1.0));
        nodes1.push_back(new Node<2>(3, false,  2.0, -1.0));

        // Create element
        std::vector<VertexElement<2,2>*> elements1;
        elements1.push_back(new VertexElement<2,2>(0, nodes1));

        // Create mesh
        VertexMesh<2,2> mesh1(nodes1, elements1);

        // Test GetShortAxisOfElement() method
        c_vector<double, 2> short_axis = mesh1.GetShortAxisOfElement(0);
        TS_ASSERT_DELTA(short_axis(0), 0.0, 1e-6);
        TS_ASSERT_DELTA(short_axis(1), 1.0, 1e-6);


        // Second test
        
        // Create nodes
        std::vector<Node<2>*> nodes2;

        // This is a rectangle, centre (0,0), width 1, length sqrt(3), rotated by 30 degrees anticlockwise                      
        nodes2.push_back(new Node<2>(0, false,  1.0, 0.0));
        nodes2.push_back(new Node<2>(1, false,  0.5, sqrt(3.0)/2.0));
        nodes2.push_back(new Node<2>(2, false, -1.0, 0.0));
        nodes2.push_back(new Node<2>(3, false, -0.5, -sqrt(3.0)/2.0));

        // Create element
        std::vector<VertexElement<2,2>*> elements2;
        elements2.push_back(new VertexElement<2,2>(0, nodes2));

        // Create mesh
        VertexMesh<2,2> mesh2(nodes2, elements2);

        // Test GetShortAxisOfElement() method
        short_axis = mesh2.GetShortAxisOfElement(0);
        TS_ASSERT_DELTA(short_axis(0), 0.5, 1e-6);
        TS_ASSERT_DELTA(short_axis(1), -sqrt(3.0)*0.5, 1e-6);


        // Third test
               
        // Test on a regular polygon (generates a random vector)
        std::vector<Node<2>*> nodes3;
        unsigned num_nodes = 6;   // vertices
        for (unsigned i=0; i<num_nodes; i++)
        {
            double theta = 2.0*M_PI*(double)(i)/(double)(num_nodes); 
            nodes3.push_back(new Node<2>(i, false, cos(theta), sin(theta)));   
        }

        // Create element
        std::vector<VertexElement<2,2>*> elements3;
        elements3.push_back(new VertexElement<2,2>(0, nodes3));

        // Create mesh
        VertexMesh<2,2> mesh3(nodes3, elements3);

        // Test GetShortAxisOfElement() method
        short_axis = mesh3.GetShortAxisOfElement(0);        
        TS_ASSERT_DELTA(short_axis(0)*short_axis(0)+short_axis(1)*short_axis(1), 1.0, 1e-6);
    }


    void TestScale(void)
    {
        // Create mesh
        VertexMesh<2,2> mesh(3, 3, 0.01, 2.0);
        
        TS_ASSERT_DELTA(mesh.GetWidth(0), 2.8867, 1e-4);
        TS_ASSERT_DELTA(mesh.GetWidth(1), 3.5000, 1e-4);

        // Squash in the x direction by a factor of 2
        mesh.Scale(0.5);

        TS_ASSERT_DELTA(mesh.GetWidth(0), 1.4433, 1e-4);
        TS_ASSERT_DELTA(mesh.GetWidth(1), 3.5000, 1e-4);

        // Stretch in the x and y directions by a factor of 2
        mesh.Scale(2.0, 2.0);

        TS_ASSERT_DELTA(mesh.GetWidth(0), 2.8867, 1e-4);
        TS_ASSERT_DELTA(mesh.GetWidth(1), 7.0000, 1e-4);
    }

};    

#endif /*TESTVERTEXMESH_HPP_*/
