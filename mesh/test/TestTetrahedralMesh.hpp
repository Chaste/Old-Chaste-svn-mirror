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

#ifndef _TESTTETRAHEDRALMESH_HPP_
#define _TESTTETRAHEDRALMESH_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>
#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "TrianglesMeshWriter.hpp"
#include "RandomNumberGenerator.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PetscTools.hpp"
#include <cmath>
#include <vector>

class TestTetrahedralMesh : public CxxTest::TestSuite
{

public:

    void TestMeshConstructionFromMeshReader(void)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 543U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 984U);

        // Check some node co-ordinates
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0],  0.9980267283, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[1], -0.0627905195, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(1)->GetPoint()[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(1)->GetPoint()[1], 0.0, 1e-6);

        // Check first element has the right nodes
        TetrahedralMesh<2,2>::ElementIterator it = mesh.GetElementIteratorBegin();
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(0), 309U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(1), 144U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(2), 310U);
        TS_ASSERT_EQUALS((*it)->GetNode(1), mesh.GetNode(144));
    }

    void TestMeshConstructionFromMeshReaderIndexedFromOne(void)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements_indexed_from_1");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 543U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 984U);

        // Check some node co-ordinates
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0], 1.0, 1e-6); // note this mesh is different to disk_984_elements
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(1)->GetPoint()[0], 0.9980267283, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(1)->GetPoint()[1], 0.0627905195, 1e-6);

        // Check first element has the right nodes
        TetrahedralMesh<2,2>::ElementIterator it = mesh.GetElementIteratorBegin();
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(0), 309U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(1), 144U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(2), 310U);
        TS_ASSERT_EQUALS((*it)->GetNode(1), mesh.GetNode(144));
    }

    void Test3dMeshConstructionFromMeshReader(void)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 51U);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 136U);
        TS_ASSERT_EQUALS(mesh_reader.GetNumFaces(), 96U);


        TetrahedralMesh<3,3> mesh;

        try
        {
            mesh.ConstructFromMeshReader(mesh_reader);
        }
        catch (Exception &e)
        {
            std::cout << e.GetMessage() << std::endl;
        }

        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 51U);
        //TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 543);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 136U);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 96U);
        TS_ASSERT_DELTA(mesh.CalculateVolume(), 1.0, 1e-15);
        TS_ASSERT_DELTA(mesh.CalculateSurfaceArea(), 6.0, 1e-16);

        // Check some node co-ordinates
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[2], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(19)->GetPoint()[0], 0.75, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(19)->GetPoint()[1], 0.25, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(19)->GetPoint()[2], 0.0, 1e-6);
    }

    void Test3dMeshConstructionFromMeshReader2(void)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_0_to_.5mm_1889_elements_irregular");
        TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 425U);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 1889U);
        TS_ASSERT_EQUALS(mesh_reader.GetNumFaces(), 436U);

        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_DELTA(mesh.CalculateVolume(), 1.25e-4, 1e-16);
        TS_ASSERT_DELTA(mesh.CalculateSurfaceArea(), 0.015, 1e-15);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 425U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 1889U);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 436U);
    }


    void TestMeshWithBoundaryElements(void)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check for the right number of boundary edges
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 100U);

        // Check all boundary elements have nodes on the boundary
        TetrahedralMesh<2,2>::BoundaryElementIterator it =
            mesh.GetBoundaryElementIteratorBegin();
        while (it != mesh.GetBoundaryElementIteratorEnd())
        {
            for (unsigned i=0; i<(*it)->GetNumNodes(); i++)
            {
                TS_ASSERT((*it)->GetNode(i)->IsBoundaryNode());
            }
            it++;
        }
    }




    void Test1DClosedMeshIn2DSpace()
    {
        TrianglesMeshReader<1,2> mesh_reader("mesh/test/data/circle_outline");
        TetrahedralMesh<1,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS( mesh.GetNumNodes(), 100U);
        TS_ASSERT_EQUALS( mesh.GetNumElements(), 100U);
        //Check that the mesh_reader has the unculled "faces" (which are nodes)
        TS_ASSERT_EQUALS(mesh_reader.GetNumFaces(), mesh.GetNumNodes());
        TS_ASSERT_EQUALS( mesh.GetNumBoundaryElements(), 0U);
    }


    void Test1DMeshIn2DSpace()
    {
        TrianglesMeshReader<1,2> mesh_reader("mesh/test/data/semicircle_outline");
        TetrahedralMesh<1,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS( mesh.GetNumNodes(), 51U);
        TS_ASSERT_EQUALS( mesh.GetNumElements(), 50U);
        //Check that the mesh_reader has the unculled "faces" (which are nodes)
        TS_ASSERT_EQUALS(mesh_reader.GetNumFaces(), mesh.GetNumNodes());
        //Culled "faces"
        TS_ASSERT_EQUALS( mesh.GetNumBoundaryElements(), 2U);
    }


    void Test2DClosedMeshIn3DSpace()
    {
        TrianglesMeshReader<2,3> mesh_reader("mesh/test/data/slab_395_elements");
        TetrahedralMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS( mesh.GetNumNodes(), 132U);
        TS_ASSERT_EQUALS( mesh.GetNumElements(), 224U);
        TS_ASSERT_EQUALS( mesh.GetNumBoundaryElements(), 0U);
    }


    void Test2DMeshIn3DSpace()
    {
        TrianglesMeshReader<2,3> mesh_reader("mesh/test/data/disk_in_3d");
        TetrahedralMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS( mesh.GetNumNodes(), 312U);
        TS_ASSERT_EQUALS( mesh.GetNumElements(), 522U);

        // Check that the mesh_reader has the culled "faces" (which are edges) (100 instead of 833)
        TS_ASSERT_EQUALS( mesh_reader.GetNumFaces(), 100U);
        // These are the 100 edges around the perimeter of the circle
        TS_ASSERT_EQUALS( mesh.GetNumBoundaryElements(), 100U);
    }


    void Test1DMeshCrossReference()
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        Node<1>* p_node = mesh.GetNode(0);
        TS_ASSERT_EQUALS(p_node->GetNumContainingElements(), 1u);
        TS_ASSERT_EQUALS(p_node->GetNumBoundaryElements(), 1u);
        Node<1>::ContainingElementIterator elt_iter = p_node->ContainingElementsBegin();
        Node<1>::ContainingBoundaryElementIterator b_elt_iter = p_node->ContainingBoundaryElementsBegin();
        TS_ASSERT_EQUALS(*elt_iter, 0u);
        TS_ASSERT_EQUALS(*b_elt_iter, 0u);

        Element<1,1>* p_element = mesh.GetElement(*elt_iter);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),0U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),1U);
        TS_ASSERT_DELTA(p_element->CalculateJacobianDeterminant(), 0.1, 1e-6);

        Node<1>* p_node2 = mesh.GetNode(1);
        TS_ASSERT_EQUALS(p_node2->GetNumContainingElements(), 2u);
        TS_ASSERT_EQUALS(p_node2->GetNumBoundaryElements(), 0u);

        elt_iter = p_node2->ContainingElementsBegin();
        p_element = mesh.GetElement(*elt_iter);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),0U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),1U);

        p_element = mesh.GetElement(*(++elt_iter));
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),1U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),2U);
        TS_ASSERT_DELTA(p_element->CalculateJacobianDeterminant(), 0.1, 1e-6);

        // There should be no more containing elements
        TS_ASSERT_EQUALS(++elt_iter, p_node2->ContainingElementsEnd());
    }


    void Test2DMeshCrossReference()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        Node<2> *p_node = mesh.GetNode(234);
        TS_ASSERT_EQUALS(p_node->GetNumContainingElements(), 5u);
        TS_ASSERT_EQUALS(p_node->GetNumBoundaryElements(), 0u);

        Node<2>::ContainingElementIterator elt_iter = p_node->ContainingElementsBegin();
        Element<2,2> *p_element;

        p_element = mesh.GetElement(*elt_iter);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),474U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),290U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2),234U);

        p_element = mesh.GetElement(*(++elt_iter));
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),234U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),461U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2),460U);

        p_element = mesh.GetElement(*(++elt_iter));
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),290U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),459U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2),234U);

        p_element = mesh.GetElement(*(++elt_iter));
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),459U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),461U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2),234U);

        p_element = mesh.GetElement(*(++elt_iter));
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),460U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),474U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2),234U);

        // Now look at a boundary node
        p_node = mesh.GetNode(99);
        TS_ASSERT_EQUALS(p_node->GetNumContainingElements(), 3u);
        TS_ASSERT_EQUALS(p_node->GetNumBoundaryElements(), 2u);
        const BoundaryElement<1,2>* p_boundary_element;
        Node<2>::ContainingBoundaryElementIterator b_elt_iter = p_node->ContainingBoundaryElementsBegin();

        p_boundary_element = mesh.GetBoundaryElement(*b_elt_iter);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(0),98U);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(1),99U);
        p_boundary_element = mesh.GetBoundaryElement(*(++b_elt_iter));
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(0),99U);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(1),0U);
    }


    void Test3DMeshCrossReference()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        Node<3>* p_node = mesh.GetNode(34);
        TS_ASSERT_EQUALS(p_node->GetNumContainingElements(), 10u);

        Element<3,3>* p_element;
        Node<3>::ContainingElementIterator elt_iter = p_node->ContainingElementsBegin();

        p_element = mesh.GetElement(*elt_iter);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),22U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),34U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2),33U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(3),10U);

        p_element = mesh.GetElement(*(++elt_iter));
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),22U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),35U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2),33U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(3),34U);

        //Now look at a boundary node
        TS_ASSERT_EQUALS(p_node->GetNumBoundaryElements(), 4u);
        const BoundaryElement<2,3> *p_boundary_element;
        Node<3>::ContainingBoundaryElementIterator b_elt_iter = p_node->ContainingBoundaryElementsBegin();
        p_boundary_element = mesh.GetBoundaryElement(*b_elt_iter);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(0),6U);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(1),34U);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(2),24U);
        p_boundary_element = mesh.GetBoundaryElement(*(++b_elt_iter));
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(0),6U);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(1),30U);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(2),34U);
        p_boundary_element = mesh.GetBoundaryElement(*(++b_elt_iter));
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(0),24U);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(1),34U);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(2),10U);
        p_boundary_element = mesh.GetBoundaryElement(*(++b_elt_iter));
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(0),34U);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(1),30U);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(2),10U);
    }

    void TestNodePermutation()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_1626_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        double volume=mesh.CalculateVolume();
        double surface=mesh.CalculateSurfaceArea();

        Node<3>*  p_node0=mesh.GetNode(0);
        Node<3>*  p_node121=mesh.GetNode(121);
        Node<3>*  p_node125=mesh.GetNode(125);
        Node<3>*  p_node273=mesh.GetNode(273);

        RandomNumberGenerator::Instance();
        mesh.PermuteNodes();

        TS_ASSERT_EQUALS(mesh.GetNode(  0)->GetIndex(),   0U);
        TS_ASSERT_EQUALS(mesh.GetNode(121)->GetIndex(), 121U);
        TS_ASSERT_EQUALS(mesh.GetNode(125)->GetIndex(), 125U);
        TS_ASSERT_EQUALS(mesh.GetNode(273)->GetIndex(), 273U);

        TS_ASSERT_EQUALS(p_node0->GetIndex(), 357U);
        TS_ASSERT_EQUALS(p_node121->GetIndex(), 35U);
        TS_ASSERT_EQUALS(p_node125->GetIndex(), 219U);
        TS_ASSERT_EQUALS(p_node273->GetIndex(), 319U);
        TS_ASSERT_EQUALS(mesh.GetNode(p_node0->GetIndex()), p_node0);
        TS_ASSERT_EQUALS(mesh.GetNode(p_node121->GetIndex()), p_node121);
        TS_ASSERT_EQUALS(mesh.GetNode(p_node125->GetIndex()), p_node125);
        TS_ASSERT_EQUALS(mesh.GetNode(p_node273->GetIndex()), p_node273);

        TS_ASSERT_DELTA(volume, mesh.CalculateVolume(), 1e-7);
        TS_ASSERT_DELTA(surface, mesh.CalculateSurfaceArea(), 1e-7);

        RandomNumberGenerator::Destroy();
    }


    void TestConstructRectangle()
    {
        TetrahedralMesh<2,2> mesh;
        unsigned width=39;
        unsigned height=16;

        mesh.ConstructRectangularMesh(width,height);

        TS_ASSERT_DELTA(mesh.CalculateVolume(), width*height, 1e-7);
        TS_ASSERT_DELTA(mesh.CalculateSurfaceArea(), 2.0*(width+height), 1e-7);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), ((width+1)*(height+1)));
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 2*(width + height));
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(),  2*(width+height) );
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2*width*height);

        TrianglesMeshWriter<2,2> mesh_writer("","RectangleMesh");
        mesh_writer.WriteFilesUsingMesh(mesh);
    }


    void TestConstructRectangleNoStagger()
    {
        TetrahedralMesh<2,2> mesh;
        unsigned width=39;
        unsigned height=16;
        mesh.ConstructRectangularMesh(width,height,false);
        TS_ASSERT_DELTA(mesh.CalculateVolume(), width*height, 1e-7);
        TS_ASSERT_DELTA(mesh.CalculateSurfaceArea(), 2.0*(width+height), 1e-7);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), ((width+1)*(height+1)));
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(),  2*(width+height) );
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2*width*height );

        TrianglesMeshWriter<2,2> mesh_writer("","RectangleMeshNoStagger");
        mesh_writer.WriteFilesUsingMesh(mesh);
    }


    void TestConstruct1x1RectangularMesh(void)
    {
        TetrahedralMesh<2,2> rect_mesh;
        rect_mesh.ConstructRectangularMesh(1, 1, false);
    }

    void TestConstructLine()
    {
        TetrahedralMesh<1,1> mesh;
        unsigned width=39;;

        mesh.ConstructLinearMesh(width);

        TS_ASSERT_DELTA(mesh.CalculateVolume(), width, 1e-7);
        TS_ASSERT_DELTA(mesh.CalculateSurfaceArea(), 0u, 1e-7);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), width+1);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 2u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(),  2u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), width);

        TrianglesMeshWriter<1,1> mesh_writer("","LineMesh");
        mesh_writer.WriteFilesUsingMesh(mesh);
    }

    void TestSetOwnerships()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        unsigned lo = 300;
        unsigned hi = 302;

        mesh.SetElementOwnerships(lo, hi);

        for (unsigned ele_num=0; ele_num< mesh.GetNumElements(); ele_num++)
        {
            bool owned = mesh.GetElement(ele_num)->GetOwnership();
            if (ele_num==26  ||
                ele_num==195 ||
                ele_num==330 ||
                ele_num==351 ||
                ele_num==498 ||
                ele_num==499 || //...these contain node 300
                ele_num==186 ||
                ele_num==208 ||
                ele_num==480 ||
                ele_num==500 ||
                ele_num==501)  //... these contain node 301
            {
                TS_ASSERT_EQUALS(owned, true);
            }
            else
            {
                TS_ASSERT_EQUALS(owned, false);
            }
        }
    }


    void TestOutwardNormal3D()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        for (unsigned i=0; i<mesh.GetNumBoundaryElements(); i++)
        {
            BoundaryElement<2,3> *b_element=mesh.GetBoundaryElement(i);
            c_vector<double, 3> normal=*(b_element->pGetWeightedDirection());
            c_vector<double, 3> centroid=b_element->CalculateCentroid();
            ChastePoint<3> out(centroid+normal);
            ChastePoint<3> in(centroid-normal);
            TS_ASSERT_THROWS_NOTHING(mesh.GetContainingElementIndex(in));
            TS_ASSERT_THROWS_ANYTHING(mesh.GetContainingElementIndex(out));
        }
    }


    void TestConstructCuboid()
    {

        TetrahedralMesh<3,3> mesh;
        unsigned width=7;
        unsigned height=4;
        unsigned depth=5;

        unsigned num_boundary_nodes =   2*( (width+1)*(height+1) + (width+1)*(depth+1) + (depth+1)*(height+1) )
                                      - 4*(width-1 + height-1 + depth-1)
                                      - 16;

        mesh.ConstructCuboid(width,height,depth);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), ((width+1)*(height+1)*(depth+1)));
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), num_boundary_nodes);

        TS_ASSERT_DELTA(mesh.CalculateVolume(), width*height*depth, 1e-7);
        TS_ASSERT_DELTA(mesh.CalculateSurfaceArea(), 2.0*(width*height+height*depth+depth*width), 1e-7);
        //Each unit square on the surface is split into 2
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(),  4*(width*height+height*depth+depth*width) );
        //Assuming that each cube is split into 6 tetrahedra
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 6*width*height*depth );

        for (unsigned i=0; i<mesh.GetNumBoundaryElements(); i++)
        {
            BoundaryElement<2,3> *b_element=mesh.GetBoundaryElement(i);
            c_vector<double, 3> normal=*(b_element->pGetWeightedDirection());
            c_vector<double, 3> centroid=b_element->CalculateCentroid();
            ChastePoint<3> out(centroid+normal);
            ChastePoint<3> in(centroid-normal);
            normal/=norm_2(normal);
            if (fabs(centroid[0]) < 1e-5)
            {
                TS_ASSERT_DELTA( normal[0], -1.0, 1e-16);
            }
            if (fabs(centroid[0] - width) < 1e-5)
            {
                TS_ASSERT_DELTA( normal[0], 1.0, 1e-16);
            }
            if (fabs(centroid[1]) < 1e-5)
            {
                TS_ASSERT_DELTA( normal[1], -1.0, 1e-16);
            }
            if (fabs(centroid[1] - height) < 1e-5)
            {
                TS_ASSERT_DELTA( normal[1], 1.0, 1e-16);
            }
            if (fabs(centroid[2]) < 1e-5)
            {
                TS_ASSERT_DELTA( normal[2], -1.0, 1e-16);
            }
            if (fabs(centroid[2] - depth) < 1e-5)
            {
                TS_ASSERT_DELTA( normal[2], 1.0, 1e-16);
            }
            TS_ASSERT_THROWS_NOTHING(mesh.GetContainingElementIndex(in));
            TS_ASSERT_THROWS_ANYTHING(mesh.GetContainingElementIndex(out));
        }

        TrianglesMeshWriter<3,3> mesh_writer("","CuboidMesh");
        mesh_writer.WriteFilesUsingMesh(mesh);
    }


    void TestPermute()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNode(0)->rGetLocation()[0], 0.0);
        TS_ASSERT_EQUALS(mesh.GetNode(0)->rGetLocation()[1], 0.0);
        TS_ASSERT_EQUALS(mesh.GetNode(0)->rGetLocation()[2], 0.0);

        TS_ASSERT_EQUALS(mesh.GetNode(1)->rGetLocation()[0], 0.2);
        TS_ASSERT_EQUALS(mesh.GetNode(1)->rGetLocation()[1], 0.0);
        TS_ASSERT_EQUALS(mesh.GetNode(1)->rGetLocation()[2], 0.0);

        TS_ASSERT_EQUALS(mesh.GetNode(2)->rGetLocation()[0], 0.2);
        TS_ASSERT_EQUALS(mesh.GetNode(2)->rGetLocation()[1], 0.2);
        TS_ASSERT_EQUALS(mesh.GetNode(2)->rGetLocation()[2], 0.0);

        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(0)->GetIndex(), 11U);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(1)->GetIndex(), 3U);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(2)->GetIndex(), 8U);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(3)->GetIndex(), 0U);

        //Make identity permuation
        std::vector<unsigned> perm;

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            perm.push_back(i);
        }
        //perm is now the identity permuation

        //Rotate first three
        perm[0]=1;
        perm[1]=2;
        perm[2]=0;

        // Rotate nodes in the first element
        perm[8]=11;
        perm[9]=8;
        perm[11]=9;

        mesh.PermuteNodes(perm);

        TS_ASSERT_EQUALS(mesh.GetNode(0)->GetIndex(), 0U);
        TS_ASSERT_EQUALS(mesh.GetNode(7)->GetIndex(), 7U);

        // 1 was node 0
        TS_ASSERT_EQUALS(mesh.GetNode(1)->rGetLocation()[0], 0.0);
        TS_ASSERT_EQUALS(mesh.GetNode(1)->rGetLocation()[1], 0.0);
        TS_ASSERT_EQUALS(mesh.GetNode(1)->rGetLocation()[2], 0.0);

        // 2 was node 1
        TS_ASSERT_EQUALS(mesh.GetNode(2)->rGetLocation()[0], 0.2);
        TS_ASSERT_EQUALS(mesh.GetNode(2)->rGetLocation()[1], 0.0);
        TS_ASSERT_EQUALS(mesh.GetNode(2)->rGetLocation()[2], 0.0);

        // 0 was node 2
        TS_ASSERT_EQUALS(mesh.GetNode(0)->rGetLocation()[0], 0.2);
        TS_ASSERT_EQUALS(mesh.GetNode(0)->rGetLocation()[1], 0.2);
        TS_ASSERT_EQUALS(mesh.GetNode(0)->rGetLocation()[2], 0.0);

        // Element 0 new indexes
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(0)->GetIndex(), 9U);  // 9 was 11
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(1)->GetIndex(), 3U);  // 3 is 3
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(2)->GetIndex(), 11U); // 11 was 8
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetNode(3)->GetIndex(), 1U);  // 1 was 0
    }


    void TestPermuteWithMetisBinaries() throw(Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_DELTA(mesh.GetNode(0)->rGetLocation()[0],  0.9980, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(0)->rGetLocation()[1], -0.0627, 1e-4);
        mesh.PermuteNodesWithMetisBinaries(4);
        TS_ASSERT_DELTA(mesh.GetNode(236)->rGetLocation()[0],  0.9980, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(236)->rGetLocation()[1], -0.0627, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(0)->rGetLocation()[0], -0.7705, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(0)->rGetLocation()[1], 0.6374, 1e-4);

        TrianglesMeshReader<3,3> mesh_reader2("mesh/test/data/3D_0_to_.5mm_1889_elements_irregular");
        TetrahedralMesh<3,3> mesh2;
        mesh2.ConstructFromMeshReader(mesh_reader2);

        TS_ASSERT_DELTA(mesh2.GetNode(0)->rGetLocation()[0], 0.0000, 1e-4);
        TS_ASSERT_DELTA(mesh2.GetNode(0)->rGetLocation()[1], 0.0000, 1e-4);
        TS_ASSERT_DELTA(mesh2.GetNode(0)->rGetLocation()[2], 0.0000, 1e-4);
        mesh2.PermuteNodesWithMetisBinaries(4);
        TS_ASSERT_DELTA(mesh2.GetNode(0)->rGetLocation()[0], 0.0000, 1e-4);
        TS_ASSERT_DELTA(mesh2.GetNode(0)->rGetLocation()[1], 0.0000, 1e-4);
        TS_ASSERT_DELTA(mesh2.GetNode(0)->rGetLocation()[2], 0.0000, 1e-4);

        TrianglesMeshWriter<3,3> mesh_writer("","3D_0_to_.5mm_1889_elements_irregular_metis");
        mesh_writer.WriteFilesUsingMesh(mesh2);
    }




    void TestClear()
    {
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(2,3);

        TS_ASSERT_EQUALS(mesh.CalculateVolume(), 6.0);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 12u);

        mesh.Clear();

        TS_ASSERT_EQUALS(mesh.GetNumNodes(),0u);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(),0u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(),0u);
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(),0u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(),0u);
        TS_ASSERT_EQUALS(mesh.GetNumAllBoundaryElements(),0u);
    }


    void TestUnflagAllElements()
    {
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(1,1);

        mesh.GetElement(0)->Flag();
        mesh.GetElement(1)->Flag();

        TS_ASSERT_EQUALS(mesh.GetElement(0)->IsFlagged(), true);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->IsFlagged(), true);

        mesh.UnflagAllElements();

        TS_ASSERT_EQUALS(mesh.GetElement(0)->IsFlagged(), false);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->IsFlagged(), false);
    }


    void TestCalculateBoundaryOfFlaggedRegion()
    {
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(3,3);

        // uncomment to write the mesh
        //TrianglesMeshWriter<2,2> mesh_writer("rectangle", "small");
        //mesh_writer.WriteFilesUsingMesh(mesh);

        mesh.GetElement(1)->Flag();
        mesh.GetElement(3)->Flag();
        mesh.GetElement(5)->Flag();
        mesh.GetElement(8)->Flag();
        mesh.GetElement(9)->Flag();
        mesh.GetElement(10)->Flag();
        mesh.GetElement(14)->Flag();
        mesh.GetElement(16)->Flag();
        mesh.GetElement(17)->Flag();

        std::set<unsigned> correct_boundary;
        correct_boundary.insert(0);
        correct_boundary.insert(4);
        correct_boundary.insert(5);
        correct_boundary.insert(2);
        correct_boundary.insert(9);
        correct_boundary.insert(13);
        correct_boundary.insert(10);
        correct_boundary.insert(7);
        correct_boundary.insert(11);
        correct_boundary.insert(14);
        correct_boundary.insert(15);

        std::set<unsigned> boundary = mesh.CalculateBoundaryOfFlaggedRegion();

        TS_ASSERT_EQUALS(correct_boundary, boundary);
    }

    void TestFlagElementsNotContainingNodes()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::set<unsigned> nodes;
        nodes.insert(0);
        nodes.insert(1);

        mesh.FlagElementsNotContainingNodes(nodes);

        TS_ASSERT_EQUALS( mesh.GetElement(0)->IsFlagged(), false);
        TS_ASSERT_EQUALS( mesh.GetElement(1)->IsFlagged(), false);
        TS_ASSERT_EQUALS( mesh.GetElement(2)->IsFlagged(), false);
        TS_ASSERT_EQUALS( mesh.GetElement(3)->IsFlagged(), true);
    }


    void TestCalculateBoundaryOfFlaggedRegion3D()
    {
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(4,4,4);
        mesh.Translate(-2,-2,-2);

        // flag elements in the positive octant
        TetrahedralMesh<3,3>::ElementIterator iter = mesh.GetElementIteratorBegin();
        while (iter != mesh.GetElementIteratorEnd())
        {
            c_vector<double, 3> centroid = (*iter)->CalculateCentroid();
            if (centroid(0)>=0 && centroid(1)>=0 && centroid(2)>=0)
            {
                (*iter)->Flag();
            }
            iter++;
        }

        // calculate boundary
        std::set<unsigned> boundary=mesh.CalculateBoundaryOfFlaggedRegion();

        // determine correct boundary
        std::set<unsigned> correct_boundary;
        for (unsigned node_index=0; node_index<mesh.GetNumNodes(); node_index++)
        {
            // node is on boundary if
            // a) 1 coordinate is zero and rest +ve or 0
            // b) 1 coordinate is 2.0 and rest +ve or 0
            // so get a sorted list of coordinates
            std::vector<double> coordinates;
            for (unsigned i=0; i<3; i++)
            {
                coordinates.push_back(mesh.GetNode(node_index)->rGetLocation()[i]);
            }
            std::sort(coordinates.begin(), coordinates.end());

            if ( (coordinates[0]==0.0)
               ||(coordinates[0]>=0.0 && coordinates[2]==2.0))
            {
                correct_boundary.insert(node_index);
            }
        }
        TS_ASSERT_EQUALS(boundary, correct_boundary);
    }


    void TestGetVectorBetweenPoints() throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        c_vector<double, 3> location1 = mesh.GetNode(0)->rGetLocation();
        c_vector<double, 3> location2 = mesh.GetNode(2)->rGetLocation();

        // test a normal distance calculation
        c_vector<double, 3> vector = mesh.GetVectorFromAtoB(location1, location2);
        TS_ASSERT_DELTA(vector[0], 0.2, 1e-7);
        TS_ASSERT_DELTA(vector[1], 0.2, 1e-7);
        TS_ASSERT_DELTA(vector[2], 0.0, 1e-7)
        TS_ASSERT_DELTA(norm_2(vector), sqrt(0.08), 1e-7);
        TS_ASSERT_DELTA(mesh.GetDistanceBetweenNodes(0, 2), sqrt(0.08), 1e-7);

        // and the opposite vector
        vector = mesh.GetVectorFromAtoB(location2, location1);
        TS_ASSERT_DELTA(vector[0], -0.2, 1e-7);
        TS_ASSERT_DELTA(vector[1], -0.2, 1e-7);
        TS_ASSERT_DELTA(vector[2], 0.0, 1e-7);
        TS_ASSERT_DELTA(norm_2(vector), sqrt(0.08), 1e-7);

        // a 3d vector
        location1[0] = 0.5;
        location1[1] = 3.0;
        location1[2] = 1.0;
        location2[0] = 2.5;
        location2[1] = 4.0;
        location2[2] = -3.0;
        vector = mesh.GetVectorFromAtoB(location1, location2);
        TS_ASSERT_DELTA(vector[0], +2.0, 1e-7);
        TS_ASSERT_DELTA(vector[1], +1.0, 1e-7);
        TS_ASSERT_DELTA(vector[2], -4.0, 1e-7);
        TS_ASSERT_DELTA(norm_2(vector), sqrt(21.0), 1e-7);
    }


    void TestMeshGetWidthAndWidthExtremesMethod(void)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");

        TetrahedralMesh<2,2> mesh;

        mesh.ConstructFromMeshReader(mesh_reader);

        double width = mesh.GetWidth(0u);
        double height = mesh.GetWidth(1u);

        TS_ASSERT_DELTA(width, 2, 1e-6);
        TS_ASSERT_DELTA(height, 2, 1e-6);

        c_vector<double,2> width_extremes = mesh.GetWidthExtremes(0u);
        c_vector<double,2> height_extremes = mesh.GetWidthExtremes(1u);

        TS_ASSERT_DELTA(width_extremes[0], -1, 1e-6);
        TS_ASSERT_DELTA(height_extremes[0], -1, 1e-6);
        TS_ASSERT_DELTA(width_extremes[1], 1, 1e-6);
        TS_ASSERT_DELTA(height_extremes[1], 1, 1e-6);
    }


    void TestPointWeightsInElement1D()
    {
        std::vector<Node<1>*> nodes1d;
        nodes1d.push_back(new Node<1>(0, false, 2.0));
        nodes1d.push_back(new Node<1>(1, false, 2.5));
        Element<1,1> element1d(INDEX_IS_NOT_USED, nodes1d);

        ChastePoint<1> in_point(2.25);
        ChastePoint<1> on_point(2.00);
        ChastePoint<1> out_point(1.25);
        c_vector <double, 2> weights;
        weights=element1d.CalculateInterpolationWeights(on_point);
        TS_ASSERT_EQUALS(weights[0], 1.0);
        TS_ASSERT_EQUALS(weights[1], 0.0);
        weights=element1d.CalculateInterpolationWeights(in_point);
        TS_ASSERT_EQUALS(weights[0], 0.5);
        TS_ASSERT_EQUALS(weights[1], 0.5);
        weights=element1d.CalculateInterpolationWeights(out_point);
        //1.25 = 2.5*2 -1.5 * 2.5
        TS_ASSERT_EQUALS(weights[0], 2.5);
        TS_ASSERT_EQUALS(weights[1], -1.5);

        delete nodes1d[0];
        delete nodes1d[1];
    }


    void TestPointInElement1D()
    {
        std::vector<Node<1>*> nodes1d;
        nodes1d.push_back(new Node<1>(0, false, 2.0));
        nodes1d.push_back(new Node<1>(1, false, 2.5));
        Element<1,1> element1d(INDEX_IS_NOT_USED, nodes1d);

        ChastePoint<1> in_point(2.25);
        ChastePoint<1> on_point(2.00);
        ChastePoint<1> out_point(1.25);
        bool strict=true;
        TS_ASSERT_EQUALS(element1d.IncludesPoint(in_point), true);
        TS_ASSERT_EQUALS(element1d.IncludesPoint(on_point), true);
        TS_ASSERT_EQUALS(element1d.IncludesPoint(on_point, strict), false);
        TS_ASSERT_EQUALS(element1d.IncludesPoint(out_point), false);

        delete nodes1d[0];
        delete nodes1d[1];
    }


    void TestPointinMesh1D(void)
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        ChastePoint<1> point1(0.15);
        ChastePoint<1> point2(-0.1);
        ChastePoint<1> point3(0.2);
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point1),1U);
        TS_ASSERT_THROWS_ANYTHING(mesh.GetContainingElementIndex(point2));
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point3),1U);  //in elements 1 and 2

        std::vector<unsigned> indices;
        indices=mesh.GetContainingElementIndices(point1);
        TS_ASSERT_EQUALS(indices.size(), 1U);
        TS_ASSERT_EQUALS(indices[0], 1U);

        indices=mesh.GetContainingElementIndices(point2);
        TS_ASSERT_EQUALS(indices.size(), 0U);

        indices=mesh.GetContainingElementIndices(point3);
        TS_ASSERT_EQUALS(indices.size(), 2U);
        TS_ASSERT_EQUALS(indices[0], 1U);
        TS_ASSERT_EQUALS(indices[1], 2U);
    }


    void TestPointWeightsAndInclusion2D()
    {
        std::vector<Node<2>*> nodes2d;
        nodes2d.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes2d.push_back(new Node<2>(1, false, 2.0, 1.0));
        nodes2d.push_back(new Node<2>(2, false, 0.0, 3.0));
        Element<2,2> element2d(INDEX_IS_NOT_USED, nodes2d);

        ChastePoint<2> on_point(0., 2.);
        c_vector <double, 3> weights;
        bool strict=true;
        TS_ASSERT_EQUALS(element2d.IncludesPoint(on_point), true);
        TS_ASSERT_EQUALS(element2d.IncludesPoint(on_point, strict), false);
        weights=element2d.CalculateInterpolationWeights(on_point);
        c_vector <double, 2> psi_on = element2d.CalculatePsi(on_point);
        TS_ASSERT_DELTA(weights[0], 1.0/3.0, 1e-5);
        TS_ASSERT_DELTA(weights[1], 0.0, 1e-5);
        TS_ASSERT_DELTA(weights[2], 2.0/3.0, 1e-5);
        TS_ASSERT_DELTA(psi_on[0],  0.0, 1e-5);
        TS_ASSERT_DELTA(psi_on[1],  2.0/3.0, 1e-5);

        ChastePoint<2> in_point(1., 1.);
        TS_ASSERT_EQUALS(element2d.IncludesPoint(in_point), true);
        weights=element2d.CalculateInterpolationWeights(in_point);
        c_vector <double, 2> psi_in = element2d.CalculatePsi(in_point);
        TS_ASSERT_LESS_THAN(0.0, weights[0]);
        TS_ASSERT_LESS_THAN(0.0, weights[1]);
        TS_ASSERT_LESS_THAN(0.0, weights[2]);
        TS_ASSERT_DELTA(psi_in[0],0.5,1e-12);
        TS_ASSERT_DELTA(psi_in[1],1.0/6.0,1e-12);

        ChastePoint<2> out_point(1., 0.);
        TS_ASSERT_EQUALS(element2d.IncludesPoint(out_point), false);
        weights=element2d.CalculateInterpolationWeights(out_point);
        c_vector <double, 2> psi_out = element2d.CalculatePsi(out_point);
        TS_ASSERT_LESS_THAN(0.0, weights[0]);
        TS_ASSERT_LESS_THAN(0.0, weights[1]);
        TS_ASSERT_LESS_THAN(weights[2], 0.0);
        TS_ASSERT_DELTA(psi_out[0],0.5,1e-12);
        TS_ASSERT_DELTA(psi_out[1],-1.0/6.0,1e-12);

        delete nodes2d[0];
        delete nodes2d[1];
        delete nodes2d[2];
    }


    void TestPointinMesh2D(void)
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_200_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        ChastePoint<2> point1(0.051, 0.051);
        ChastePoint<2> point2(0.2,0.2);
        ChastePoint<2> point3(0.05, 0.05); //Node 60 of mesh
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point1),110U);
        TS_ASSERT_EQUALS(mesh.GetNearestElementIndex(point1),110U);
        TS_ASSERT_THROWS_ANYTHING(mesh.GetContainingElementIndex(point2));
        TS_ASSERT_EQUALS(mesh.GetNearestElementIndex(point2),199U); //Contains top-right corner
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point3),89U);  //in elements 89,90,91,108,109, 110

        std::vector<unsigned> indices;
        indices=mesh.GetContainingElementIndices(point1);
        TS_ASSERT_EQUALS(indices.size(), 1U);
        TS_ASSERT_EQUALS(indices[0], 110U);

        indices=mesh.GetContainingElementIndices(point2);
        TS_ASSERT_EQUALS(indices.size(), 0U);

        indices=mesh.GetContainingElementIndices(point3);
        TS_ASSERT_EQUALS(indices.size(), 6U);
        TS_ASSERT_EQUALS(indices[0], 89U);
        TS_ASSERT_EQUALS(indices[1], 90U);
        TS_ASSERT_EQUALS(indices[5], 110U);
    }


    void TestPointInElement3D()
    {
        std::vector<Node<3>*> nodes3d;
        nodes3d.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        nodes3d.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        nodes3d.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        nodes3d.push_back(new Node<3>(3, false, 0.0, 0.0, 1.0));
        Element<3,3> element3d(INDEX_IS_NOT_USED, nodes3d);

        bool strict=true;
        ChastePoint<3> on_point(0., 0.2, 0.);
        c_vector <double, 4> weights;
        TS_ASSERT_EQUALS(element3d.IncludesPoint(on_point), true);
        TS_ASSERT_EQUALS(element3d.IncludesPoint(on_point, strict), false);
        weights=element3d.CalculateInterpolationWeights(on_point);
        c_vector <double, 3> psi_on = element3d.CalculatePsi(on_point);
        TS_ASSERT_DELTA(weights[0], 0.8, 1e-5);
        TS_ASSERT_DELTA(weights[1], 0.0, 1e-5);
        TS_ASSERT_DELTA(weights[2], 0.2, 1e-5);
        TS_ASSERT_DELTA(weights[3], 0.0, 1e-5);
        TS_ASSERT_DELTA(psi_on[0],0.0,1e-12);
        TS_ASSERT_DELTA(psi_on[1],0.2,1e-12);
        TS_ASSERT_DELTA(psi_on[2],0.0,1e-12);

        ChastePoint<3> in_point(.25, .25, .25);
        TS_ASSERT_EQUALS(element3d.IncludesPoint(in_point), true);
        weights=element3d.CalculateInterpolationWeights(in_point);
        TS_ASSERT_LESS_THAN(0.0, weights[0]);
        TS_ASSERT_LESS_THAN(0.0, weights[1]);
        TS_ASSERT_LESS_THAN(0.0, weights[2]);
        TS_ASSERT_LESS_THAN(0.0, weights[3]);
        // Weights are non-negative and sum to 1
        TS_ASSERT_DELTA(norm_1(weights),1.0,1e-12);

        ChastePoint<3> out_point(0.1, -10., 0.1);
        TS_ASSERT_EQUALS(element3d.IncludesPoint(out_point), false);
        weights=element3d.CalculateInterpolationWeights(out_point);
        TS_ASSERT_LESS_THAN(0.0, weights[0]);
        TS_ASSERT_LESS_THAN(0.0, weights[1]);
        TS_ASSERT_LESS_THAN(weights[2], 0.0);
        TS_ASSERT_LESS_THAN(0.0, weights[3]);
        TS_ASSERT_DELTA(weights[0], 10.8, 1e-5);
        TS_ASSERT_DELTA(weights[1], 0.1, 1e-5);
        TS_ASSERT_DELTA(weights[2], -10.0, 1e-5);
        TS_ASSERT_DELTA(weights[3], 0.1, 1e-5);
        // Weights still sum to 1, but one weight is negative
        TS_ASSERT_DELTA(norm_1(weights),21.0,1e-12);

        weights=element3d.CalculateInterpolationWeightsWithProjection(out_point);
        TS_ASSERT_LESS_THAN(0.0, weights[0]);
        TS_ASSERT_LESS_THAN(0.0, weights[1]);
        TS_ASSERT_LESS_THAN_EQUALS(0.0, weights[2]);
        TS_ASSERT_EQUALS(0.0, weights[2]);
        TS_ASSERT_LESS_THAN(0.0, weights[3]);
        // Weights are non-negative and sum to 1
        TS_ASSERT_DELTA(norm_1(weights),1.0,1e-12);

        delete nodes3d[0];
        delete nodes3d[1];
        delete nodes3d[2];
        delete nodes3d[3];
    }


    void TestPointinMesh3D(void)
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_0_to_1mm_6000_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        ChastePoint<3> point1(0.051, 0.051,0.051);
        ChastePoint<3> point2(0.2,0.2,0.2);
        ChastePoint<3> point3(0.050000000000000003,  0.050000000000000003,  0.050000000000000003);
        //Node 665 of mesh
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point1),2992U);
        TS_ASSERT_THROWS_ANYTHING(mesh.GetContainingElementIndex(point2));
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point3),2044U);
        /*in elements 2044, 2047. 2058, 2192, 2268, 2286, 2392, 2414, 2415,
         * 2424, 2426, 2452, 2661, 2704, 2734, 2745, 2846, 2968, 2990, 2992,
         * 3015, 3022, 3024, 3026
         */

        //should throw because vertex is not strictly contained in any element
        TS_ASSERT_THROWS_ANYTHING(mesh.GetContainingElementIndex(point3, true));

        std::vector<unsigned> indices;
        indices=mesh.GetContainingElementIndices(point1);
        TS_ASSERT_EQUALS(indices.size(), 1U);
        TS_ASSERT_EQUALS(indices[0], 2992U);

        indices=mesh.GetContainingElementIndices(point2);
        TS_ASSERT_EQUALS(indices.size(), 0U);

        indices=mesh.GetContainingElementIndices(point3);
        TS_ASSERT_EQUALS(indices.size(), 24U);
        TS_ASSERT_EQUALS(indices[0], 2044U);
        TS_ASSERT_EQUALS(indices[1], 2047U);
        TS_ASSERT_EQUALS(indices[5], 2286U);
        TS_ASSERT_EQUALS(indices[23], 3026U);

        // test when a suggested set of elements is given
        std::set<unsigned> suggested_elements;
        suggested_elements.insert(2991);
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point1,false,suggested_elements),2992U);
        suggested_elements.insert(2992); // should find it quicker, but can't really test that
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point1,false,suggested_elements),2992U);
        suggested_elements.insert(2993);
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point1,false,suggested_elements),2992U);
    }


    void TestFloatingPointIn3D()
    {
        //There's some weird failing behaviour in the refined mesh test
        //This test duplicates it

        TetrahedralMesh<3,3> mesh;

        mesh.ConstructCuboid(3, 3, 3);
        double third=1.0L/3.0L;
        mesh.Scale(third, third, third);

        ChastePoint<3> point_on_edge1(5.0/6.0,   0.5,       1.0);
        ChastePoint<3> point_on_edge2(5.0L/6.0L, 0.5,       1.0);
        ChastePoint<3> point_on_edge3(5.0L/6.0L, 0.5L,      1.0L);
        ChastePoint<3> point_on_edge4(5.0L/6.0L, 3.0L/6.0L, 1.0L);
        ChastePoint<3> point_on_edge5(5.0L/6.0L, 0.5L,      6.0L/6.0L);
        ChastePoint<3> point_on_edge6(5.0L/6.0L, 3.0L/6.0L, 6.0L/6.0L);
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point_on_edge1),142U);
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point_on_edge2),142U);
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point_on_edge3),142U);
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point_on_edge4),142U);
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point_on_edge5),142U);
        TS_ASSERT_EQUALS(mesh.GetContainingElementIndex(point_on_edge6),142U);
    }


    /*
     * This tests that a 'dummy' archive function does not throw any errors
     * (a mesh writer stores the mesh in a nice format anyway, we only
     * need this so that subclasses can archive their own member variables
     */
    void TestArchiveTetrahedralMesh()
    {
        EXIT_IF_PARALLEL;
        OutputFileHandler handler("archive", false);    // do not clean folder
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "conf_mesh.arch";

        // Create an ouput archive
        {
            TetrahedralMesh<2,2>* const p_mesh = new TetrahedralMesh<2,2>;
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
            TetrahedralMesh<2,2> mesh;
            p_mesh->ConstructFromMeshReader(mesh_reader);
            TS_ASSERT_EQUALS(p_mesh->GetNumNodes(),4u);
            TS_ASSERT_EQUALS(p_mesh->GetNumElements(),2u);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_mesh;
            delete p_mesh;
        }

        {
            TetrahedralMesh<2,2>* p_mesh2;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // restore from the archive
            input_arch >> p_mesh2;

            // \TODO:ticket:412 These lines will test proper mesh archiving.
//            TS_ASSERT_EQUALS(p_mesh2->GetNumNodes(),4u);
//            TS_ASSERT_EQUALS(p_mesh2->GetNumElements(),2u);

            delete p_mesh2;
        }
    }


    void TestGetAngleBetweenNodes()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_DELTA( mesh.GetAngleBetweenNodes(0,1),  0.0,      1e-12);
        TS_ASSERT_DELTA( mesh.GetAngleBetweenNodes(0,2),  M_PI/4,   1e-12);
        TS_ASSERT_DELTA( mesh.GetAngleBetweenNodes(0,3),  M_PI/2,   1e-12);
        TS_ASSERT_DELTA( mesh.GetAngleBetweenNodes(1,0),  M_PI,     1e-12);
        TS_ASSERT_DELTA( mesh.GetAngleBetweenNodes(2,0), -3*M_PI/4, 1e-12);

        TS_ASSERT_THROWS_ANYTHING(mesh.GetAngleBetweenNodes(0,0));
    }


    void TestNodesPerProcessorFile() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // throws because file does not exist
        TS_ASSERT_THROWS_ANYTHING(mesh.ReadNodesPerProcessorFile("dsgund"));

        // throws because sum of nodes is not equal to the number of nodes in the mesh
        TS_ASSERT_THROWS_ANYTHING(mesh.ReadNodesPerProcessorFile("mesh/test/data/nodes_per_processor_1.txt"));

        mesh.ReadNodesPerProcessorFile("mesh/test/data/nodes_per_processor_2.txt");

        // nodes_per_processor_2.txt = {1,3}
        TS_ASSERT_EQUALS(mesh.rGetNodesPerProcessor().size(), 2u);
        TS_ASSERT_EQUALS(mesh.rGetNodesPerProcessor()[0], 1u);
        TS_ASSERT_EQUALS(mesh.rGetNodesPerProcessor()[1], 3u);
    }

    void TestReadingMeshesWithRegions() throw (Exception)
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements_with_attributes");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 1U);

        for (unsigned i=0; i<10; i++)
        {
            TS_ASSERT_EQUALS(mesh.GetElement(i)->GetRegion(), i%5+1);
        }
    }
};
#endif //_TESTTETRAHEDRALMESH_HPP_
