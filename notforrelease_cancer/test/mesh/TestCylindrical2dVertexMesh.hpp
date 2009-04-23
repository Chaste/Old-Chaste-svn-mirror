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
#ifndef TESTCYLINDRICAL2DVERTEXMESH_HPP_
#define TESTCYLINDRICAL2DVERTEXMESH_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "Cylindrical2dVertexMesh.hpp"
#include "VertexMeshWriter.hpp"

class TestCylindrical2dVertexMesh : public CxxTest::TestSuite
{
public:

    void TestMeshGenerator()
    {
        // Create non-periodic mesh
        VertexMesh<2,2> non_periodic_mesh(4, 4, 0.01, 2.0);

        TS_ASSERT_EQUALS(non_periodic_mesh.GetNumElements(), 16u);
        TS_ASSERT_EQUALS(non_periodic_mesh.GetNumNodes(), 48u);

        // Create periodic mesh
        Cylindrical2dVertexMesh cylindrical_vertex_mesh(4, 4, 0.01, 2.0);

        // The periodic mesh should have the same number of elements but fewer nodes
        TS_ASSERT_EQUALS(cylindrical_vertex_mesh.GetNumElements(), 16u);
        TS_ASSERT_EQUALS(cylindrical_vertex_mesh.GetNumNodes(), 40u);

        // Create a vertex mesh writer with cylindrical mesh
        VertexMeshWriter<2,2> vertex_mesh_writer("TestCylindrical2dVertexMesh", "cylindrical_vertex_mesh");
        vertex_mesh_writer.WriteFilesUsingMesh(cylindrical_vertex_mesh);

        OutputFileHandler handler("TestCylindrical2dVertexMesh", false);
        std::string results_file1 = handler.GetOutputDirectoryFullPath() + "cylindrical_vertex_mesh.node";
        std::string results_file2 = handler.GetOutputDirectoryFullPath() + "cylindrical_vertex_mesh.cell";

        TS_ASSERT_EQUALS(system(("diff " + results_file1 + " notforrelease_cancer/test/data/TestCylindrical2dVertexMesh/cylindrical_vertex_mesh.node").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_file2 + " notforrelease_cancer/test/data/TestCylindrical2dVertexMesh/cylindrical_vertex_mesh.cell").c_str()), 0);
    }


    void TestMeshGetWidth()
    {
        // Create mesh
        Cylindrical2dVertexMesh cylindrical_vertex_mesh(4, 4, 0.01, 2.0);

        // Test GetWidthExtremes() method
        c_vector<double,2> width_extremes = cylindrical_vertex_mesh.GetWidthExtremes(0u);
        c_vector<double,2> height_extremes = cylindrical_vertex_mesh.GetWidthExtremes(1u);

        TS_ASSERT_DELTA(width_extremes[0], 0.0000, 1e-4);
        TS_ASSERT_DELTA(width_extremes[1], 2.8867, 1e-4);

        TS_ASSERT_DELTA(height_extremes[0], 0.0000, 1e-4);
        TS_ASSERT_DELTA(height_extremes[1], 4.5000, 1e-4);

        // Test GetWidth() method
        double width = cylindrical_vertex_mesh.GetWidth(0);
        double height = cylindrical_vertex_mesh.GetWidth(1);

        TS_ASSERT_DELTA(width, 3.4641, 1e-4);
        TS_ASSERT_DELTA(height, 4.5000, 1e-4);
    }


    void TestGetVectorFromAtoB() throw (Exception)
    {
        // Create mesh
        Cylindrical2dVertexMesh mesh(4, 4, 0.01, 2.0);

        c_vector<double, 2> node18_location = mesh.GetNode(18)->rGetLocation();
        c_vector<double, 2> node19_location = mesh.GetNode(19)->rGetLocation();

        // Test a normal vector and distance calculation
        c_vector<double, 2> vector = mesh.GetVectorFromAtoB(node18_location, node19_location);
        TS_ASSERT_DELTA(vector[0], 0.5773, 1e-4);
        TS_ASSERT_DELTA(vector[1], 0.0000, 1e-4);
        TS_ASSERT_DELTA(norm_2(vector), 0.5773, 1e-4);
        TS_ASSERT_DELTA(mesh.GetDistanceBetweenNodes(18, 19), 0.5773, 1e-4);

        // Test the opposite vector
        vector = mesh.GetVectorFromAtoB(node19_location, node18_location);
        TS_ASSERT_DELTA(vector[0], -0.5773, 1e-4);
        TS_ASSERT_DELTA(vector[1], 0.0000, 1e-4);

        // Test a periodic calculation
        c_vector<double, 2> node16_location = mesh.GetNode(16)->rGetLocation();
        vector = mesh.GetVectorFromAtoB(node16_location, node19_location);

        TS_ASSERT_DELTA(vector[0], -1.1547, 1e-4);
        TS_ASSERT_DELTA(vector[1], 0.0000, 1e-4);
    }


    void TestSetNodeLocationForCylindricalMesh() throw (Exception)
    {
        // Create mesh
        Cylindrical2dVertexMesh mesh(4, 4, 0.01, 2.0);

        // Move one of the nodes to near the periodic boundary
        c_vector<double, 2> new_point_location;
        new_point_location[0] = -0.01;
        new_point_location[1] = 1.5;
        ChastePoint<2> new_point(new_point_location);

        // This node was on left and is now near the right
        mesh.SetNode(12u, new_point);
        TS_ASSERT_DELTA(mesh.GetNode(12u)->rGetLocation()[0], 3.4541, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(12u)->rGetLocation()[1], 1.5, 1e-4);

        // This node has stayed close to where it was
        new_point.SetCoordinate(0u, 0.2);
        mesh.SetNode(0u, new_point);
        TS_ASSERT_DELTA(mesh.GetNode(0u)->rGetLocation()[0], 0.2, 1e-4);

        // This node was on right and is now near the left
        new_point.SetCoordinate(0u, 3.5);
        mesh.SetNode(8u, new_point);
        TS_ASSERT_DELTA(mesh.GetNode(8u)->rGetLocation()[0], 0.0358, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(8u)->rGetLocation()[1], 1.5, 1e-4);
    }


    void TestAddNodeAndReMesh() throw (Exception)
    {
        // Create mesh
        Cylindrical2dVertexMesh mesh(6, 6, 0.01, 2.0);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 84u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 36u);

        // Choose a node on the left boundary
        ChastePoint<2> point = mesh.GetNode(18)->GetPoint();
        TS_ASSERT_DELTA(point[0], 0.0000, 1e-4);
        TS_ASSERT_DELTA(point[1], 1.5000, 1e-4);

        // Create a new node close to this node
        point.SetCoordinate(0, -0.01);
        point.SetCoordinate(1, 4.5);
        Node<2>* p_node = new Node<2>(mesh.GetNumNodes(), point);

        unsigned old_num_nodes = mesh.GetNumNodes();

        // Add this new node to the mesh
        unsigned new_index = mesh.AddNode(p_node);
        TS_ASSERT_EQUALS(new_index, old_num_nodes);

        // Remesh to update correspondences
        mesh.SetEdgeDivisionThreshold(1000); // set high threshold to avoid more nodes appearing in the remesh
        VertexElementMap map(mesh.GetNumElements());
        mesh.ReMesh(map);

        TS_ASSERT_EQUALS(map.Size(), mesh.GetNumElements());
        TS_ASSERT_EQUALS(map.IsIdentityMap(), true);

        // Check that the mesh is updated
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 85u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 36u);

        TS_ASSERT_DELTA(mesh.GetNode(new_index)->rGetLocation()[0], 5.18615, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(new_index)->rGetLocation()[1], 4.5000, 1e-4);

        // Now tet AddNode() when mDeletedNodeIndices is populated

        // Label node 29 as deleted
        mesh.mDeletedNodeIndices.push_back(29);

        // Create a new node close to this node
        ChastePoint<2> point2;
        point2.SetCoordinate(0, 2.0);
        point2.SetCoordinate(1, 2.1);
        Node<2>* p_node2 = new Node<2>(mesh.GetNumNodes(), point);

        // Add this new node to the mesh
        new_index = mesh.AddNode(p_node2);
        TS_ASSERT_EQUALS(new_index, 29u);
    }


    void TestElementAreaPerimeterAndCentroid()
    {
        // Create mesh
        Cylindrical2dVertexMesh mesh(4, 4, 0.01, 2.0);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 40u);

        // Test area and perimeter calculations for  all  elements
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            TS_ASSERT_DELTA(mesh.GetAreaOfElement(i), 0.8660, 1e-4);
            TS_ASSERT_DELTA(mesh.GetPerimeterOfElement(i), 3.4641, 1e-4);
        }

        // Test centroid calculations for nonperiodic element
        c_vector<double, 2> centroid = mesh.GetCentroidOfElement(5);
        TS_ASSERT_DELTA(centroid(0), 1.4433, 1e-4);
        TS_ASSERT_DELTA(centroid(1), 2.0, 1e-4);

        // Test centroid calculations for periodic element
        centroid = mesh.GetCentroidOfElement(7);
        TS_ASSERT_DELTA(centroid(0), 3.1754, 1e-4);
        TS_ASSERT_DELTA(centroid(1), 2.0, 1e-4);
    }


    void TestArchiving() throw (Exception)
    {
        std::string dirname = "archive";
        OutputFileHandler handler(dirname, false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "cylindrical_vertex_mesh_base.arch";

        std::string mesh_filename = "cylindrical_vertex_mesh";
        std::string mesh_pathname = handler.GetOutputDirectoryFullPath() + mesh_filename;

        // Create mesh
        unsigned num_cells_across = 4;
        unsigned num_cells_up = 7;

        Cylindrical2dVertexMesh* const p_mesh = new Cylindrical2dVertexMesh(num_cells_across, num_cells_up, 0.01, 2.0);

        double crypt_width = num_cells_across*sqrt(3.0)/2.0;

        /*
         * You need the const above to stop a BOOST_STATIC_ASSERTION failure.
         * This is because the serialization library only allows you to save
         * tracked objects while the compiler considers them const, to prevent
         * the objects changing during the save, and so object tracking leading
         * to wrong results. For example, A is saved once via pointer, then
         * changed, then saved again.  The second save notes that A was saved
         * before, so doesn't write its data again, and the change is lost.
         */
        {
            // Serialize the mesh
            double width = p_mesh->GetWidth(0);
            TS_ASSERT_DELTA(width, crypt_width, 1e-7);

            // Save the mesh data using mesh writers
            VertexMeshWriter<2,2> vertex_mesh_writer(dirname, mesh_filename);
            vertex_mesh_writer.WriteFilesUsingMesh(*p_mesh);

            // Archive the mesh
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // We have to serialize via a pointer here, or the derived class information is lost.
            output_arch << p_mesh;
        }

        {
            // De-serialize and compare
            Cylindrical2dVertexMesh* p_mesh2;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_mesh2;

            // Re-initialise the mesh
            p_mesh2->Clear();
            VertexMeshReader<2,2> mesh_reader(mesh_pathname);
            p_mesh2->ConstructFromMeshReader(mesh_reader);

            // Compare the loaded mesh against the original

            // Compare width
            TS_ASSERT_DELTA(p_mesh2->GetWidth(0), crypt_width, 1e-7);

            // Compare nodes
            TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), p_mesh2->GetNumNodes());

            for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
            {
                Node<2>* p_node = p_mesh->GetNode(i);
                Node<2>* p_node2 = p_mesh2->GetNode(i);
                TS_ASSERT_EQUALS(p_node->IsDeleted(), p_node2->IsDeleted());
                TS_ASSERT_EQUALS(p_node->GetIndex(), p_node2->GetIndex());
                TS_ASSERT_EQUALS(p_node->IsBoundaryNode(), p_node2->IsBoundaryNode());
                for (unsigned j=0; j<2; j++)
                {
                    TS_ASSERT_DELTA(p_node->rGetLocation()[j], p_node2->rGetLocation()[j], 1e-4);
                }
            }

            // Compare elements
            TS_ASSERT_EQUALS(p_mesh->GetNumElements(), p_mesh2->GetNumElements());
            TS_ASSERT_EQUALS(p_mesh->GetNumAllElements(), p_mesh2->GetNumAllElements());

            for (unsigned i=0; i<p_mesh->GetNumElements(); i++)
            {
                VertexElement<2,2>* p_elt = p_mesh->GetElement(i);
                VertexElement<2,2>* p_elt2 = p_mesh2->GetElement(i);
                TS_ASSERT_EQUALS(p_elt->GetNumNodes(), p_elt2->GetNumNodes());
                for (unsigned i=0; i<p_elt->GetNumNodes(); i++)
                {
                    TS_ASSERT_EQUALS(p_elt->GetNodeGlobalIndex(i), p_elt2->GetNodeGlobalIndex(i));
                }
            }

            // Tidy up
            delete p_mesh;
            delete p_mesh2;
        }
    }


    void TestCylindricalReMesh() throw (Exception)
    {
        // Create mesh
        unsigned num_cells_across = 6;
        unsigned num_cells_up = 12;

        Cylindrical2dVertexMesh mesh(num_cells_across, num_cells_up, 0.01, 2.0);

        // Remesh
        VertexElementMap map(mesh.GetNumElements());
        mesh.ReMesh(map);

        TS_ASSERT_EQUALS(map.Size(), mesh.GetNumElements());
        TS_ASSERT_EQUALS(map.IsIdentityMap(), true);

        // Check that there are the correct number of everything
        TS_ASSERT_EQUALS(mesh.GetNumElements(), num_cells_across*num_cells_up);
    }


    void TestCylindricalReMeshAfterDelete() throw (Exception)
    {
        // Create mesh
        unsigned num_cells_across = 6;
        unsigned num_cells_up = 12;

        Cylindrical2dVertexMesh mesh(num_cells_across, num_cells_up, 0.01, 2.0);
        unsigned num_old_nodes = mesh.GetNumNodes();
        unsigned num_old_elements = num_cells_across*num_cells_up;

        // Delete a node
        mesh.DeleteElementPriorToReMesh(8);

        // Remesh
        VertexElementMap map(mesh.GetNumElements());
        mesh.ReMesh(map);

        TS_ASSERT_EQUALS(map.IsIdentityMap(), false);
        TS_ASSERT_EQUALS(map.Size(), num_old_elements);

        // Check that there are the correct number of everything
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), num_old_nodes);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), num_old_elements-1);
   }

};

#endif /*TESTCYLINDRICAL2DVERTEXMESH_HPP_*/
