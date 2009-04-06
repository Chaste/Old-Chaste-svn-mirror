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


#ifndef TESTREMESH_HPP_
#define TESTREMESH_HPP_

#include <cxxtest/TestSuite.h>
#include "MutableMesh.hpp"
#include <cmath>


class TestRemesh : public CxxTest::TestSuite
{
public:

    /**
     * Test 3D remesh - very similar test to TestOperationOfTetgenMoveNodes,
     * but uses mesh.Remesh() instead of calling tetgen from here.
     */
    void TestRemesh3dMoveNodes() throw (Exception)
    {
        OutputFileHandler handler("");

        TrianglesMeshReader<3,3> mesh_reader2("mesh/test/data/cube_1626_elements");
        TetrahedralMesh<3,3> old_mesh;
        old_mesh.ConstructFromMeshReader(mesh_reader2);

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_1626_elements");
        MutableMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            ChastePoint<3> point = mesh.GetNode(i)->GetPoint();
            ChastePoint<3> old_mesh_point = old_mesh.GetNode(i)->GetPoint();

            for (int j=0; j<3; j++)
            {
                if (fabs(point[j]-0.0) >1e-6 && fabs(point[j]-1.0) >1e-6)
                {
                    point.SetCoordinate(j, point[j]+9e-2);
                    old_mesh_point.SetCoordinate(j, old_mesh_point[j]+9e-2);
                }
            }

            mesh.GetNode(i)->SetPoint(point);
            old_mesh.GetNode(i)->SetPoint(old_mesh_point);
        }
        mesh.RefreshMesh();
        old_mesh.RefreshMesh();

        double old_volume = mesh.CalculateVolume();
        TS_ASSERT_DELTA(1, old_volume, 1e-7);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 375u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 1626u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 390u);

        NodeMap map(mesh.GetNumNodes());
        mesh.ReMesh(map);

        TS_ASSERT_EQUALS(map.Size(), mesh.GetNumNodes());

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), old_mesh.GetNumNodes());
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), old_mesh.GetNumBoundaryElements());

        TS_ASSERT_EQUALS(mesh.GetNumElements()+1, old_mesh.GetNumElements());

        // Test to see whether triangle/ tetgen is renumbering the nodes
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            // The map turns out to be the identity map in this test
            TS_ASSERT_EQUALS(map.GetNewIndex(i), i);

            const c_vector<double, 3> node_loc1 = mesh.GetNode(map.GetNewIndex(i))->rGetLocation();
            const c_vector<double, 3> node_loc2 = old_mesh.GetNode(i)->rGetLocation();

            for (int j=0; j<3; j++)
            {
                TS_ASSERT_DELTA(node_loc1[j], node_loc2[j], 1e-6);
            }
        }

        double new_volume = mesh.CalculateVolume();
        TS_ASSERT_DELTA(old_volume, new_volume, 1e-7);
    }

    void TestOperationOfTetgenMoveNodes() throw (Exception)
    {
        OutputFileHandler handler("");

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_1626_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            ChastePoint<3> point = mesh.GetNode(i)->GetPoint();

            for (int j=0; j<3; j++)
            {
                if (fabs(point[j]-0.0) >1e-6 && fabs(point[j]-1.0) >1e-6)
                {
                    point.SetCoordinate(j, point[j]+9e-2);
                }
            }

            mesh.GetNode(i)->SetPoint(point);
        }
        mesh.RefreshMesh();

        double volume = mesh.CalculateVolume();
        TS_ASSERT_DELTA(1, volume, 1e-7);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 375u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 1626u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 390u);

        out_stream node_file = handler.OpenOutputFile("temp.node");
        (*node_file) << mesh.GetNumNodes() << "\t3\t0\t0\n";
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            const c_vector<double, 3> node_loc = mesh.GetNode(i)->rGetLocation();
            (*node_file) << i << "\t" << node_loc[0] << "\t" << node_loc[1] << "\t" << node_loc[2] << "\n";
        }
        node_file->close();
        std::string full_name = handler.GetOutputDirectoryFullPath("")+"temp.";
        std::string command = "tetgen -Q " + full_name + "node" + " > /dev/null";
        int return_value=system(command.c_str());
        TS_ASSERT(return_value == 0);
        if (return_value != 0)
        {
            TS_TRACE("Do you have tetgen from http://tetgen.berlios.de/ installed in your path?");
        }
        else
        {
            TrianglesMeshReader<3,3> mesh_reader2(full_name+"1");
            TetrahedralMesh<3,3> mesh2;
            mesh2.ConstructFromMeshReader(mesh_reader2);
            TS_ASSERT_EQUALS(mesh.GetNumNodes(), mesh2.GetNumNodes());
            TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), mesh2.GetNumBoundaryElements());

            TS_ASSERT_EQUALS(mesh.GetNumElements(), mesh2.GetNumElements()+1);

            // Test to see whether triangle/tetgen is renumbering the nodes

            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                const c_vector<double, 3> node_loc1 = mesh.GetNode(i)->rGetLocation();
                const c_vector<double, 3> node_loc2 = mesh2.GetNode(i)->rGetLocation();

                for (int j=0; j<3; j++)
                {
                    TS_ASSERT_DELTA(node_loc1[j], node_loc2[j], 1e-6);
                }
            }
        }
    }


    void TestRemeshWithMethod2D() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        double area = mesh.CalculateVolume();
        const int node_index = 432;
        const int target_index = 206;

        unsigned num_nodes_before = mesh.GetNumNodes();
        unsigned num_elements_before = mesh.GetNumElements();
        unsigned num_boundary_elements_before = mesh.GetNumBoundaryElements();

        mesh.MoveMergeNode(node_index, target_index);

        TS_ASSERT_DELTA(area, mesh.CalculateVolume(), 1e-6);
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), mesh.GetNumElements()+2);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), mesh.GetNumNodes()+1);
        TS_ASSERT_EQUALS(mesh.GetNumAllBoundaryElements(), mesh.GetNumBoundaryElements());

        NodeMap map(1);
        mesh.ReMesh(map);

        TS_ASSERT_EQUALS(map.Size(), mesh.GetNumNodes()+1);//one node removed during remesh

        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), mesh.GetNumElements());
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(),mesh.GetNumNodes());
        TS_ASSERT_EQUALS(mesh.GetNumAllBoundaryElements(), mesh.GetNumBoundaryElements());

        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), num_elements_before-2);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), num_nodes_before-1);
        TS_ASSERT_EQUALS(mesh.GetNumAllBoundaryElements(), num_boundary_elements_before);
        TS_ASSERT_DELTA(mesh.CalculateVolume(), area, 1e-6);
    }

    void TestRemeshWithMethod3D() throw (Exception)
    {
        // Create mutable tetrahedral mesh which is Delaunay
        std::vector<Node<3> *> nodes;

        nodes.push_back(new Node<3>(0, true,  0.0,  0.0,  0.0));
        nodes.push_back(new Node<3>(1, true,  1.0,  1.0,  0.0));
        nodes.push_back(new Node<3>(2, true,  1.0,  0.0,  1.0));
        nodes.push_back(new Node<3>(3, true,  0.0,  1.0,  1.0));
        nodes.push_back(new Node<3>(4, false, 0.5,  0.5,  0.5));

        MutableMesh<3,3> mesh(nodes);
        double area = mesh.CalculateVolume();

        unsigned num_nodes_before = mesh.GetNumNodes();
        unsigned num_elements_before = mesh.GetNumElements();
        unsigned num_boundary_elements_before = mesh.GetNumBoundaryElements();

        mesh.DeleteNode(4);
        TS_ASSERT_DELTA(area, mesh.CalculateVolume(), 1e-6);

        NodeMap map(1);
        mesh.ReMesh(map);

        TS_ASSERT_EQUALS(map.Size(), mesh.GetNumNodes()+1);//one node removed during remesh

        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), mesh.GetNumElements());
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(),mesh.GetNumNodes());
        TS_ASSERT_EQUALS(mesh.GetNumAllBoundaryElements(), mesh.GetNumBoundaryElements());

        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), num_elements_before-3);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), num_nodes_before-1);
        TS_ASSERT_EQUALS(mesh.GetNumAllBoundaryElements(), num_boundary_elements_before);
        TS_ASSERT_DELTA(mesh.CalculateVolume(), area, 1e-6);

        mesh.DeleteNode(3);
        TS_ASSERT_THROWS_ANYTHING(mesh.ReMesh(map));
    }

    void TestNodeMap()
    {
        NodeMap map(10);
        TS_ASSERT_EQUALS(map.Size(), 10u);

        map.ResetToIdentity();
        TS_ASSERT_EQUALS(map.IsIdentityMap(), true);

        map.SetNewIndex(0,1);
        map.SetNewIndex(1,0);

        TS_ASSERT_EQUALS(map.GetNewIndex(0), 1u);
        TS_ASSERT_EQUALS(map.GetNewIndex(1), 0u);
        TS_ASSERT_EQUALS(map.GetNewIndex(2), 2u);

        TS_ASSERT_EQUALS(map.IsIdentityMap(), false);

        map.ResetToIdentity();
        map.SetDeleted(4);
        TS_ASSERT_EQUALS(map.IsDeleted(4), true);
        TS_ASSERT_EQUALS(map.IsDeleted(5), false);
        TS_ASSERT_EQUALS(map.IsIdentityMap(), false);
    }

    void TestReMeshFailsAfterEnoughDeletions() throw (Exception)
    {
        MutableMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(1,1);
        NodeMap map(1);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 4u);
        mesh.ReMesh(map);

        mesh.DeleteNodePriorToReMesh(3);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 4u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 3u);

        mesh.ReMesh(map);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 3u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 3u);

        mesh.DeleteNodePriorToReMesh(2);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 3u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 2u);

        TS_ASSERT_THROWS_ANYTHING(mesh.ReMesh(map));
    }

    void TestRawTriangleLibraryCall()
    {
        struct triangulateio in, out;

        /* Define input points. */

        in.numberofpoints = 5;
        in.numberofpointattributes = 0;
        in.pointlist = (double *) malloc(in.numberofpoints * 2 * sizeof(double));
        in.pointlist[0] = 0.0;
        in.pointlist[1] = 0.0;
        in.pointlist[2] = 1.0;
        in.pointlist[3] = 0.0;
        in.pointlist[4] = 1.0;
        in.pointlist[5] = 10.0;
        in.pointlist[6] = 0.0;
        in.pointlist[7] = 10.0;
        in.pointlist[8] = 0.5;
        in.pointlist[9] = 7.0;

        in.pointmarkerlist = NULL;
        in.numberofsegments = 0;
        in.numberofholes = 0;
        in.numberofregions = 0;

        out.pointlist = NULL;
        out.pointattributelist = (double *) NULL;
        out.pointmarkerlist = (int *) NULL;
        out.trianglelist = (int *) NULL;
        out.triangleattributelist = (double *) NULL;
        out.edgelist = (int *) NULL;
        out.edgemarkerlist = (int *) NULL;

        triangulate((char*)"Qze", &in, &out, NULL);

        TS_ASSERT_EQUALS(out.numberofpoints,             5);
        TS_ASSERT_EQUALS(out.numberofpointattributes,    0);
        TS_ASSERT_EQUALS(out.numberoftriangles,          4);
        TS_ASSERT_EQUALS(out.numberofcorners,            3);
        TS_ASSERT_EQUALS(out.numberoftriangleattributes, 0);
        TS_ASSERT_EQUALS(out.numberofedges,              8);

        // Free all allocated arrays, including those allocated by Triangle
        free(in.pointlist);

        free(out.pointlist);
        free(out.pointattributelist);
        free(out.pointmarkerlist);
        free(out.trianglelist);
        free(out.triangleattributelist);
        free(out.edgelist);
        free(out.edgemarkerlist);
    }

    void TestRemeshWithLibraryMethodSimple() throw (Exception)
    {
        // Same data as previous test
        std::vector<Node<2> *> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 1.0, 10.0));
        nodes.push_back(new Node<2>(3, true, 0.0, 10.0));
        nodes.push_back(new Node<2>(4, true, 0.5, 7.0));

        MutableMesh<2,2> mesh(nodes);

        TS_ASSERT_DELTA(mesh.CalculateVolume(), 10.0, 1e-6);
        TS_ASSERT_DELTA(mesh.CalculateSurfaceArea(), 22.0, 1e-6);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 4u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 4u);

        NodeMap map(1);
        mesh.ReMesh(map);

        TS_ASSERT_EQUALS(map.Size(), mesh.GetNumNodes());

        TS_ASSERT_DELTA(mesh.CalculateVolume(), 10.0, 1e-6);
        TS_ASSERT_DELTA(mesh.CalculateSurfaceArea(), 22.0, 1e-6);

        TS_ASSERT_EQUALS(mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 4u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 4u);
    }

    void TestRemeshWithLibraryMethod2D() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        double area = mesh.CalculateVolume();
        const int node_index = 432;
        const int target_index = 206;

        unsigned num_nodes_before = mesh.GetNumNodes();
        unsigned num_elements_before = mesh.GetNumElements();
        unsigned num_boundary_elements_before = mesh.GetNumBoundaryElements();

        mesh.MoveMergeNode(node_index, target_index);

        TS_ASSERT_DELTA(area, mesh.CalculateVolume(), 1e-6);
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), mesh.GetNumElements()+2);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), mesh.GetNumNodes()+1);
        TS_ASSERT_EQUALS(mesh.GetNumAllBoundaryElements(), mesh.GetNumBoundaryElements());

        NodeMap map(1);
        mesh.ReMesh(map);

        TS_ASSERT_EQUALS(map.Size(), mesh.GetNumNodes()+1); //one node removed during remesh

        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), mesh.GetNumElements());
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), mesh.GetNumNodes());
        TS_ASSERT_EQUALS(mesh.GetNumAllBoundaryElements(), mesh.GetNumBoundaryElements());

        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), num_elements_before-2);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), num_nodes_before-1);
        TS_ASSERT_EQUALS(mesh.GetNumAllBoundaryElements(), num_boundary_elements_before);
        TS_ASSERT_DELTA(mesh.CalculateVolume(), area, 1e-6);
    }

};

#endif /*TESTREMESH_HPP_*/
