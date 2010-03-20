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

        //Make the faces
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

           //Make the elements
           std::vector<VertexElement<2,3>*> faces_element_0, faces_element_1;
        std::vector<bool> orientations_0, orientations_1;

        //Cube element
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

        //Pyramid element
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

        TS_ASSERT_EQUALS(mesh.GetNumNodes(),9u);
        TS_ASSERT_EQUALS(mesh.GetNumFaces(),10u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(),2u);

        //test Location of random node
        TS_ASSERT_DELTA(mesh.GetNode(2)->rGetLocation()[0], 0.0, 1e-3);
        TS_ASSERT_DELTA(mesh.GetNode(2)->rGetLocation()[1], 1.0, 1e-3);
        TS_ASSERT_DELTA(mesh.GetNode(2)->rGetLocation()[2], 0.0, 1e-3);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(2)->GetIndex(),2u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(2)->GetIndex(),6u);

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


    void TestGetCentroidOfElement() throw(Exception)
    {
//        // Create nodes
//        std::vector<Node<2>*> nodes;
//        unsigned num_nodes = 6;
//        for (unsigned i=0; i<num_nodes; i++)
//        {
//            double theta = 2.0*M_PI*(double)(i)/(double)(num_nodes);
//            nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
//        }
//
//        // Create element
//        std::vector<VertexElement<2,2>*> elements;
//        elements.push_back(new VertexElement<2,2>(0, nodes));
//
//        // Create mesh
//        MutableVertexMesh<2,2> mesh(nodes, elements);
//
//        // Test GetCentroidOfElement() method
//        c_vector<double, 2> centroid = mesh.GetCentroidOfElement(0);
//
//        TS_ASSERT_DELTA(centroid(0), 0.0, 1e-6);
//        TS_ASSERT_DELTA(centroid(1), 0.0, 1e-6);
    }

    void TestVertexElementAreaAndPerimeter()
    {
//        // Create nodes
//        std::vector<Node<2>*> nodes;
//        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
//        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
//        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
//        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
//
//        // Create element
//        std::vector<VertexElement<2,2>*> elements;
//        elements.push_back(new VertexElement<2,2>(0, nodes));
//
//        // Create mesh
//        MutableVertexMesh<2,2> mesh(nodes, elements);
//
//        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 4u);
//
//        // Check nodes have correct indices
//        for (unsigned i=0; i<4; i++)
//        {
//            TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNodeGlobalIndex(i), i);
//        }
//
//        // Test area and perimeter calculations
//        TS_ASSERT_DELTA(mesh.GetAreaOfElement(0), 1.0, 1e-6);
//        TS_ASSERT_DELTA(mesh.GetPerimeterOfElement(0), 4.0, 1e-6);
    }


    void TestMeshGetWidthAndWidthExtremes()
    {
//        // Create mesh
//        HoneycombVertexMeshGenerator generator(3, 3);
//        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
//
//        // Test GetWidthExtremes() method
//        c_vector<double,2> width_extremes = p_mesh->GetWidthExtremes(0u);
//        c_vector<double,2> height_extremes = p_mesh->GetWidthExtremes(1u);
//
//        TS_ASSERT_DELTA(width_extremes[0], 0.0000, 1e-4);
//        TS_ASSERT_DELTA(width_extremes[1], 3.5000, 1e-4);
//
//        TS_ASSERT_DELTA(height_extremes[0], 0.0000, 1e-4);
//        TS_ASSERT_DELTA(height_extremes[1], 2.8867, 1e-4);
//
//        // Test GetWidth() method
//        double width = p_mesh->GetWidth(0);
//        double height = p_mesh->GetWidth(1);
//
//        TS_ASSERT_DELTA(height, 2.8867, 1e-4);
//        TS_ASSERT_DELTA(width, 3.5000, 1e-4);
    }

    void TestMeshConstructionFromMeshReader()
    {
//        // Create mesh
//        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/vertex_mesh");
//        MutableVertexMesh<2,2> mesh;
//        mesh.ConstructFromMeshReader(mesh_reader);
//
//        // Test Get methods
//        TS_ASSERT_DELTA(mesh.GetCellRearrangementThreshold(), 0.01, 1e-4); // Default value
//        TS_ASSERT_DELTA(mesh.GetEdgeDivisionThreshold(), DBL_MAX, 1e-4); // Default value
//        TS_ASSERT_DELTA(mesh.GetT2Threshold(), 0.001, 1e-4); // Default value
//
//        // Check we have the right number of nodes & elements
//        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 7u);
//        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);
//
//        // Check some node co-ordinates
//        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0], 0.0, 1e-6);
//        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[1], 0.0, 1e-6);
//        TS_ASSERT_DELTA(mesh.GetNode(2)->GetPoint()[0], 1.5, 1e-6);
//        TS_ASSERT_DELTA(mesh.GetNode(2)->GetPoint()[1], 1.0, 1e-6);
//
//        // Check second element has the right nodes
//        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(0), 2u);
//        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(1), 5u);
//        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(2), 6u);
//        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(1), mesh.GetNode(5));
//
//        // Create mesh in which elements have attributes
//        VertexMeshReader<2,2> mesh_reader2("mesh/test/data/TestVertexMesh/vertex_mesh_with_attributes");
//        MutableVertexMesh<2,2> mesh2;
//        mesh2.ConstructFromMeshReader(mesh_reader2);
//
//        // Check we have the right number of nodes & elements
//        TS_ASSERT_EQUALS(mesh2.GetNumNodes(), 7u);
//        TS_ASSERT_EQUALS(mesh2.GetNumElements(), 2u);
//
//        // Check some node co-ordinates
//        TS_ASSERT_DELTA(mesh2.GetNode(0)->GetPoint()[0], 0.0, 1e-6);
//        TS_ASSERT_DELTA(mesh2.GetNode(0)->GetPoint()[1], 0.0, 1e-6);
//        TS_ASSERT_DELTA(mesh2.GetNode(2)->GetPoint()[0], 1.5, 1e-6);
//        TS_ASSERT_DELTA(mesh2.GetNode(2)->GetPoint()[1], 1.0, 1e-6);
//
//        // Check second element has the right nodes
//        TS_ASSERT_EQUALS(mesh2.GetElement(1)->GetNodeGlobalIndex(0), 2u);
//        TS_ASSERT_EQUALS(mesh2.GetElement(1)->GetNodeGlobalIndex(1), 5u);
//        TS_ASSERT_EQUALS(mesh2.GetElement(1)->GetNodeGlobalIndex(2), 6u);
//        TS_ASSERT_EQUALS(mesh2.GetElement(1)->GetNode(1), mesh2.GetNode(5));
//
//        // Check element attributes
//        TS_ASSERT_EQUALS(mesh2.GetElement(0)->GetRegion(), 76u);
//        TS_ASSERT_EQUALS(mesh2.GetElement(1)->GetRegion(), 89u);
    }

    void TestMeshConstructionFromMeshReaderIndexedFromOne()
    {
//        // Create mesh
//        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMesh/vertex_mesh_elements_indexed_from_1");
//        MutableVertexMesh<2,2> mesh;
//        mesh.ConstructFromMeshReader(mesh_reader);
//
//        // Check we have the right number of nodes & elements
//        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 7u);
//        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);
//
//        // Check some node co-ordinates
//        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0], 0.0, 1e-6);
//        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[1], 0.0, 1e-6);
//        TS_ASSERT_DELTA(mesh.GetNode(2)->GetPoint()[0], 1.5, 1e-6);
//        TS_ASSERT_DELTA(mesh.GetNode(2)->GetPoint()[1], 1.0, 1e-6);
//
//        // Check first element has the right nodes
//        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(0), 2u);
//        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(1), 5u);
//        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(2), 6u);
//        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(1), mesh.GetNode(5));
    }

    // This tests that a 'dummy' archive function does not throw any errors
    void TestArchiveVertexMesh()
    {
//        std::string archive_dir = "archive";
//        std::string archive_file = "vertex_mesh_base.arch";
//        ArchiveLocationInfo::SetMeshFilename("vertex_mesh");
//
//        HoneycombVertexMeshGenerator generator(5, 3);
//        AbstractMesh<2,2>* const p_mesh = generator.GetMesh();
//
//        /*
//         * You need the const above to stop a BOOST_STATIC_ASSERTION failure.
//         * This is because the serialization library only allows you to save tracked
//         * objects while the compiler considers them const, to prevent the objects
//         * changing during the save, and so object tracking leading to wrong results.
//         *
//         * E.g. A is saved once via pointer, then changed, then saved again. The second
//         * save notes that A was saved before, so doesn't write its data again, and the
//         * change is lost.
//         */
//
//        // Create an output archive
//        {
//            TS_ASSERT_EQUALS( (static_cast<VertexMesh<2,2>* >(p_mesh))->GetNumNodes(), 46u);
//            TS_ASSERT_EQUALS( (static_cast<VertexMesh<2,2>* >(p_mesh))->GetNumElements(), 15u);
//
//            // Create output archive
//            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
//            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();
//
//            // We have to serialize via a pointer here, or the derived class information is lost
//            (*p_arch) << p_mesh;
//        }
//
//        {
//            // De-serialize and compare
//            AbstractMesh<2,2>* p_mesh2;
//
//            // Create an input archive
//            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
//            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();
//
//            // Restore from the archive
//            (*p_arch) >> p_mesh2;
//
//            MutableVertexMesh<2,2>* p_mesh_original = static_cast<MutableVertexMesh<2,2>*>(p_mesh2);
//            MutableVertexMesh<2,2>* p_mesh_loaded = static_cast<MutableVertexMesh<2,2>*>(p_mesh);
//
//            // Compare the loaded mesh against the original
//
//            TS_ASSERT_EQUALS(p_mesh_original->GetNumNodes(), p_mesh_loaded->GetNumNodes());
//
//            for (unsigned node_index=0; node_index<p_mesh_original->GetNumNodes(); node_index++)
//            {
//                Node<2>* p_node = p_mesh_original->GetNode(node_index);
//                Node<2>* p_node2 = p_mesh_loaded->GetNode(node_index);
//
//                TS_ASSERT_EQUALS(p_node->IsDeleted(), p_node2->IsDeleted());
//                TS_ASSERT_EQUALS(p_node->GetIndex(), p_node2->GetIndex());
//
/////\todo This line was commented as part of #1076 - will reinstate once reading/writing of boundary elements
/////      is done properly for vertex meshes
////                TS_ASSERT_EQUALS(p_node->IsBoundaryNode(), p_node2->IsBoundaryNode());
//
//                for (unsigned dimension=0; dimension<2; dimension++)
//                {
//                    TS_ASSERT_DELTA(p_node->rGetLocation()[dimension], p_node2->rGetLocation()[dimension], 1e-4);
//                }
//            }
//
//            TS_ASSERT_EQUALS(p_mesh_original->GetNumElements(), p_mesh_loaded->GetNumElements());
//
//            for (unsigned elem_index=0; elem_index < p_mesh_original->GetNumElements(); elem_index++)
//            {
//                TS_ASSERT_EQUALS(p_mesh_original->GetElement(elem_index)->GetNumNodes(),
//                                 p_mesh_loaded->GetElement(elem_index)->GetNumNodes());
//
//                for (unsigned local_index=0; local_index<p_mesh_original->GetElement(elem_index)->GetNumNodes(); local_index++)
//                {
//                    TS_ASSERT_EQUALS(p_mesh_original->GetElement(elem_index)->GetNodeGlobalIndex(local_index),
//                                     p_mesh_loaded->GetElement(elem_index)->GetNodeGlobalIndex(local_index));
//                }
//            }
//
//            // Tidy up
//            delete p_mesh_original;
//        }
    }

    void TestScaleAndTranslate()
    {
//        // Create 2D mesh
//        HoneycombVertexMeshGenerator generator(3, 3);
//        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
//
//        TS_ASSERT_DELTA(p_mesh->GetWidth(0), 3.5000, 1e-4);
//        TS_ASSERT_DELTA(p_mesh->GetWidth(1), 2.8867, 1e-4);
//
//        // Squash in the x direction by a factor of 2
//        p_mesh->Scale(0.5);
//
//        TS_ASSERT_DELTA(p_mesh->GetWidth(0), 1.7500, 1e-4);
//        TS_ASSERT_DELTA(p_mesh->GetWidth(1), 2.8867, 1e-4);
//
//        // Stretch in the x and y directions by a factor of 2
//        p_mesh->Scale(2.0, 2.0);
//
//        TS_ASSERT_DELTA(p_mesh->GetWidth(0), 3.5000, 1e-4);
//        TS_ASSERT_DELTA(p_mesh->GetWidth(1), 5.7735, 1e-4);
//
//        // Create 3D mesh
//        std::vector<Node<3>*> nodes;
//        nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
//        nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
//        nodes.push_back(new Node<3>(2, false, 1.0, 2.0, 0.0));
//        nodes.push_back(new Node<3>(3, false, 0.0, 2.0, 0.0));
//        nodes.push_back(new Node<3>(4, false, 0.0, 2.0, 3.0));
//        nodes.push_back(new Node<3>(5, false, 1.0, 0.0, 3.0));
//        nodes.push_back(new Node<3>(6, false, 1.0, 2.0, 3.0));
//        nodes.push_back(new Node<3>(7, false, 0.0, 2.0, 3.0));
//
//        std::vector<VertexElement<3,3>*> elements;
//        elements.push_back(new VertexElement<3,3>(0, nodes));
//
//        MutableVertexMesh<3,3> mesh3d(nodes, elements);
//
//        TS_ASSERT_DELTA(mesh3d.GetWidth(0), 1.0, 1e-4);
//        TS_ASSERT_DELTA(mesh3d.GetWidth(1), 2.0, 1e-4);
//        TS_ASSERT_DELTA(mesh3d.GetWidth(2), 3.0, 1e-4);
//
//        // Stretch the mesh
//        mesh3d.Scale(4.0, 2.0, 4.0/3.0);
//
//        TS_ASSERT_DELTA(mesh3d.GetWidth(0), 4.0, 1e-4);
//        TS_ASSERT_DELTA(mesh3d.GetWidth(1), 4.0, 1e-4);
//        TS_ASSERT_DELTA(mesh3d.GetWidth(2), 4.0, 1e-4);
//
//        // Test the translate method
//        // Pick a certain node and store spatial position
//        Node<3>* p_node = mesh3d.GetNode(7);
//        ChastePoint<3> original_coordinate = p_node->GetPoint();
//
//        const double x_movement = 1.0;
//        const double y_movement = 2.5;
//        const double z_movement = 2.5;
//
//        mesh3d.Translate(x_movement, y_movement, z_movement);
//
//        ChastePoint<3>  new_coordinate = p_node->GetPoint();
//
//        TS_ASSERT_DELTA(original_coordinate[0], new_coordinate[0] - x_movement, 1e-6);
//        TS_ASSERT_DELTA(original_coordinate[1], new_coordinate[1] - y_movement, 1e-6);
//        TS_ASSERT_DELTA(original_coordinate[2], new_coordinate[2] - z_movement, 1e-6);
    }

    void TestBoundaryNodes()
    {
//        // Create a mesh with just boundary nodes
//        std::vector<Node<2>*> nodes;
//        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
//        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
//        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
//        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
//
//        std::vector<VertexElement<2,2>*> elements;
//        elements.push_back(new VertexElement<2,2>(0, nodes));
//
//        MutableVertexMesh<2,2> mesh1(nodes, elements);
//
//        // Test boundary property of nodes
//        for (unsigned i=0; i<mesh1.GetNumNodes(); i++)
//        {
//            TS_ASSERT_EQUALS(mesh1.GetNode(i)->IsBoundaryNode(), false);
//        }
//
//        // Create a mesh with some interior nodes
//        HoneycombVertexMeshGenerator generator1(2, 2, false, 0.01, 2.0);
//        MutableVertexMesh<2,2>* p_mesh1 = generator1.GetMesh();
//
//        // Test boundary property of nodes
//        for (unsigned i=0; i<p_mesh1->GetNumNodes(); i++)
//        {
//            bool expected_boundary_node = (i==6 || i==9) ? false : true;
//            TS_ASSERT_EQUALS(p_mesh1->GetNode(i)->IsBoundaryNode(), expected_boundary_node);
//        }
//
//        // Create a larger mesh with some interior nodes
//        HoneycombVertexMeshGenerator generator2(3, 3, false, 0.01, 2.0);
//        MutableVertexMesh<2,2>* p_mesh2 = generator2.GetMesh();
//
//        // Test boundary property of nodes
//        for (unsigned i=0; i<p_mesh2->GetNumNodes(); i++)
//        {
//            bool expected_boundary_node = true;
//            if (i==8 || i==9 || i==12 || i==13 || i==16 || i==17 || i==20 || i==21)
//            {
//                expected_boundary_node = false;
//            }
//
//            TS_ASSERT_EQUALS(p_mesh2->GetNode(i)->IsBoundaryNode(), expected_boundary_node);
//        }
    }

    void TestTranslation2DWithUblas()
    {
//        // Create 2D mesh
//        HoneycombVertexMeshGenerator generator(3, 3, false, 0.01, 2.0);
//        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
//
//        c_vector<double, 2> old_location1 = p_mesh->GetNode(4)->rGetLocation();
//        c_vector<double, 2> old_location2 = p_mesh->GetNode(9)->rGetLocation();
//
//        // Set translation vector
//        c_vector<double, 2> trans_vec;
//        trans_vec(0) = 2.0;
//        trans_vec(1) = 3.0;
//
//        // Translate
//        p_mesh->Translate(trans_vec);
//        c_vector<double, 2> new_location1 = p_mesh->GetNode(4)->rGetLocation();
//        c_vector<double, 2> new_location2 = p_mesh->GetNode(9)->rGetLocation();
//
//        // Spot check a couple of nodes
//        TS_ASSERT_DELTA(new_location1[0], old_location1[0] + 2.0, 1e-6);
//        TS_ASSERT_DELTA(new_location1[1], old_location1[1] + 3.0, 1e-6);
//
//        TS_ASSERT_DELTA(new_location2[0], old_location2[0] + 2.0, 1e-6);
//        TS_ASSERT_DELTA(new_location2[1], old_location2[1] + 3.0, 1e-6);
    }

    void TestTranslation2DMethod() throw (Exception)
    {
//        // Create 2D mesh
//        HoneycombVertexMeshGenerator generator(3, 3);
//        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
//
//        // Pick a random node and store spatial position
//        Node<2>* p_node = p_mesh->GetNode(10);
//        ChastePoint<2> original_coordinate = p_node->GetPoint();
//
//        const double x_movement = 1.0;
//        const double y_movement = 2.5;
//
//        p_mesh->Translate(x_movement, y_movement);
//
//        ChastePoint<2>  new_coordinate = p_node->GetPoint();
//
//        TS_ASSERT_DELTA(original_coordinate[0], new_coordinate[0] - x_movement, 1e-6);
//        TS_ASSERT_DELTA(original_coordinate[1], new_coordinate[1] - y_movement, 1e-6);
    }

};

#endif /*TESTVERTEXMESH_HPP_*/
