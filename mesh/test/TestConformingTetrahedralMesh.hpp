#ifndef _TESTCONFORMINGTETRAHEDRALMESH_HPP_
#define _TESTCONFORMINGTETRAHEDRALMESH_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include "TrianglesMeshWriter.cpp"
#include "RandomNumberGenerator.hpp"
#include <cmath>

#include <vector>

class TestConformingTetrahedralMesh : public CxxTest::TestSuite
{
public:

    void TestMeshConstructionFromMeshReader(void)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        
        ConformingTetrahedralMesh<2,2> mesh;
        
        mesh.ConstructFromMeshReader(mesh_reader,1);
        
        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumCornerNodes(), 543);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 984);
        
        // Check some node co-ordinates
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0],  0.9980267283, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[1], -0.0627905195, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(1)->GetPoint()[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(1)->GetPoint()[1], 0.0, 1e-6);
        
        // Check first element has the right nodes
        ConformingTetrahedralMesh<2,2>::ElementIterator it = mesh.GetElementIteratorBegin();
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(0), 309);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(1), 144);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(2), 310);
        TS_ASSERT_EQUALS((*it)->GetNode(1), mesh.GetNode(144));
        
    }
    
    void TestSimpleQuadraticMeshConstructionFromMeshReader(void)
    {
    
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        
        ConformingTetrahedralMesh<2,2> mesh;
        
        try
        {
            mesh.ConstructFromMeshReader(mesh_reader,2);
        }
        catch (Exception &e)
        {
            std::cout << e.GetMessage() << std::endl;
            TS_ASSERT(0);
        }
        
        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumCornerNodes(), 4);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 9);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2);
        
        // Check some node co-ordinates
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(1)->GetPoint()[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(1)->GetPoint()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(2)->GetPoint()[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(2)->GetPoint()[1], 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(3)->GetPoint()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(3)->GetPoint()[1], 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(4)->GetPoint()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(4)->GetPoint()[1], 0.5, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(5)->GetPoint()[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(5)->GetPoint()[1], 0.5, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(6)->GetPoint()[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(6)->GetPoint()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(7)->GetPoint()[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(7)->GetPoint()[1], 0.5, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(8)->GetPoint()[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(8)->GetPoint()[1], 1.0, 1e-6);
        
        // Check all elements have the right nodes
        ConformingTetrahedralMesh<2,2>::ElementIterator it = mesh.GetElementIteratorBegin();
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(0), 3);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(1), 0);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(2), 1);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(3), 4);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(4), 5);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(5), 6);
//		TS_ASSERT_EQUALS((*it)->GetNode(1), mesh.GetNode(144));
        it++;
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(0), 1);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(1), 2);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(2), 3);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(3), 7);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(4), 5);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(5), 8);
        
    }
    
    void TestQuadraticMeshConstructionFromMeshReader(void)
    {
    
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        
        ConformingTetrahedralMesh<2,2> mesh;
        
        try
        {
            mesh.ConstructFromMeshReader(mesh_reader,2);
        }
        catch (Exception &e)
        {
            std::cout << e.GetMessage() << std::endl;
            TS_ASSERT(0);
        }
        
        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumCornerNodes(), 543);
        //TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 543);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 984);
        
        // Check some node co-ordinates
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0],  0.9980267283, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[1], -0.0627905195, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(1)->GetPoint()[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(1)->GetPoint()[1], 0.0, 1e-6);
        
        // Check first element has the right nodes
        ConformingTetrahedralMesh<2,2>::ElementIterator it = mesh.GetElementIteratorBegin();
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(0), 309);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(1), 144);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(2), 310);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(3), 543);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(4), 544);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(5), 545);
        TS_ASSERT_EQUALS((*it)->GetNode(1), mesh.GetNode(144));
        it++;
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(3), 546);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(4), 547);
        
    }
    
    void Test3dLinearMeshConstructionFromMeshReader(void)
    {
    
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        
        //const int DIM = pMeshReader->GetDimension();
        ConformingTetrahedralMesh<3,3> mesh;
        
        try
        {
            mesh.ConstructFromMeshReader(mesh_reader,1);
        }
        catch (Exception &e)
        {
            std::cout << e.GetMessage() << std::endl;
        }
        
        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumCornerNodes(), 51);
        //TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 543);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 136);
        
        // Check some node co-ordinates
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[2], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(19)->GetPoint()[0], 0.75, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(19)->GetPoint()[1], 0.25, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(19)->GetPoint()[2], 0.0, 1e-6);
        
    }
    
    void Test3dQuadraticMeshConstructionFromMeshReader(void)
    {
    
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        
        //const int DIM = pMeshReader->GetDimension();
        ConformingTetrahedralMesh<3,3> mesh;
        
        try
        {
            mesh.ConstructFromMeshReader(mesh_reader,2);
        }
        catch (Exception &e)
        {
            std::cout << e.GetMessage() << std::endl;
        }
        
        // Check we have the right number of nodes & elements
        TS_ASSERT_EQUALS(mesh.GetNumCornerNodes(), 51);
        //TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 543);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 136);
        
        // Check some node co-ordinates
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[2], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(19)->GetPoint()[0], 0.75, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(19)->GetPoint()[1], 0.25, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(19)->GetPoint()[2], 0.0, 1e-6);
        
        // Check first element has the right nodes
        ConformingTetrahedralMesh<3,3>::ElementIterator it = mesh.GetElementIteratorBegin();
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(0), 17);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(1), 10);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(2), 16);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(3), 18);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(4), 51);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(5), 52);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(6), 53);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(7), 54);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(8), 55);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(9), 56);
        TS_ASSERT_EQUALS((*it)->GetNode(5), mesh.GetNode(52));
        it++;
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(5), 58);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(6), 59);
        
    }
    
    void TestMeshWithBoundaryElements(void)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Check for the right number of boundary edges
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 100);
        
        // Check all boundary elements have nodes on the boundary
        ConformingTetrahedralMesh<2,2>::BoundaryElementIterator it =
            mesh.GetBoundaryElementIteratorBegin();
        while (it != mesh.GetBoundaryElementIteratorEnd())
        {
            for (int i=0; i<(*it)->GetNumNodes(); i++)
            {
                TS_ASSERT((*it)->GetNode(i)->IsBoundaryNode());
            }
            it++;
        }
    }
    
    void TestRescaleMeshFromBoundaryNode(void)
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        Point<1> updatedPoint(1.5);
        mesh.RescaleMeshFromBoundaryNode(updatedPoint,10);
        for (int i=0; i < 11; i++)
        {
            TS_ASSERT_DELTA(mesh.GetNode(i)->GetPoint()[0], 1.5*(i/10.0) , 0.001);
        }
    }
    
    void Test1DClosedMeshIn2DSpace()
    {
        TrianglesMeshReader<1,2> mesh_reader("mesh/test/data/circle_outline");
        ;
        ConformingTetrahedralMesh<1,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TS_ASSERT_EQUALS( mesh.GetNumNodes(), 100);
        TS_ASSERT_EQUALS( mesh.GetNumElements(), 100);
        //Check that the mesh_reader has the unculled "faces" (which are nodes)
        TS_ASSERT_EQUALS( mesh_reader.GetNumFaces(), mesh.GetNumNodes());
        TS_ASSERT_EQUALS( mesh.GetNumBoundaryElements(), 0);
        
    }
    
    void Test1DMeshIn2DSpace()
    {
        TrianglesMeshReader<1,2> mesh_reader("mesh/test/data/semicircle_outline");
        ;
        ConformingTetrahedralMesh<1,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TS_ASSERT_EQUALS( mesh.GetNumNodes(), 51);
        TS_ASSERT_EQUALS( mesh.GetNumElements(), 50);
        //Check that the mesh_reader has the unculled "faces" (which are nodes)
        TS_ASSERT_EQUALS( mesh_reader.GetNumFaces(), mesh.GetNumNodes());
        //Culled "faces"
        TS_ASSERT_EQUALS( mesh.GetNumBoundaryElements(), 2);
    }
    
    
    void Test2DClosedMeshIn3DSpace()
    {
        TrianglesMeshReader<2,3> mesh_reader("mesh/test/data/slab_395_elements");
        ConformingTetrahedralMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TS_ASSERT_EQUALS( mesh.GetNumNodes(), 132);
        TS_ASSERT_EQUALS( mesh.GetNumElements(), 224);
        TS_ASSERT_EQUALS( mesh.GetNumBoundaryElements(), 0);
    }
    
    void Test2DMeshIn3DSpace()
    {
        TrianglesMeshReader<2,3> mesh_reader("mesh/test/data/disk_in_3d");
        ConformingTetrahedralMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TS_ASSERT_EQUALS( mesh.GetNumNodes(), 312);
        TS_ASSERT_EQUALS( mesh.GetNumElements(), 522);
        
        
        ///Check that the mesh_reader has the unculled "faces" (which are edges)
        TS_ASSERT_EQUALS( mesh_reader.GetNumFaces(), 833);
        //These are the 100 edges around the perimeter of the circle
        TS_ASSERT_EQUALS( mesh.GetNumBoundaryElements(), 100);
    }
    
    
    void Test1DMeshCrossReference()
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        
        ConformingTetrahedralMesh<1,1> mesh;
        
        mesh.ConstructFromMeshReader(mesh_reader);
        
        Node<1> *p_node=mesh.GetNode(0);
        TS_ASSERT_EQUALS(p_node->GetNumContainingElements(), 1u);
        TS_ASSERT_EQUALS(p_node->GetNumBoundaryElements(), 1u);
        unsigned boundary_element_index= p_node->GetNextBoundaryElementIndex();
        TS_ASSERT_EQUALS(boundary_element_index, 0u);
        unsigned element_index= p_node->GetNextContainingElementIndex();
        TS_ASSERT_EQUALS(element_index, 0u);
        
        Element<1,1> * p_element = mesh.GetElement(element_index);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),0);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),1);
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.1, 1e-6);
        
        Node<1> *p_node2=mesh.GetNode(1);
        TS_ASSERT_EQUALS(p_node2->GetNumContainingElements(), 2u);
        TS_ASSERT_EQUALS(p_node2->GetNumBoundaryElements(), 0u);
        
        p_element = mesh.GetElement(p_node2->GetNextContainingElementIndex());
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),0);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),1);
        
        p_element =mesh.GetElement(p_node2->GetNextContainingElementIndex());
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),1);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),2);
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.1, 1e-6);
        
        // This should wrap back to 1st element
        p_element = mesh.GetElement(p_node2->GetNextContainingElementIndex());
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),0);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),1);
    }
    void Test2DMeshCrossReference()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        
        ConformingTetrahedralMesh<2,2> mesh;
        
        mesh.ConstructFromMeshReader(mesh_reader);
        
        Node<2> *p_node=mesh.GetNode(234);
        TS_ASSERT_EQUALS(p_node->GetNumContainingElements(), 5u);
        TS_ASSERT_EQUALS(p_node->GetNumBoundaryElements(), 0u);
        Element<2,2> *p_element;
        
        p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),474);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),290);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2),234);
        
        p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),234);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),461);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2),460);
        
        p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),290);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),459);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2),234);
        
        p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),459);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),461);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2),234);
        
        p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),460);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),474);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2),234);
        
        //Now look at a boundary node
        p_node=mesh.GetNode(99);
        TS_ASSERT_EQUALS(p_node->GetNumContainingElements(), 3u);
        TS_ASSERT_EQUALS(p_node->GetNumBoundaryElements(), 2u);
        const BoundaryElement<1,2> *p_boundary_element;
        p_boundary_element = mesh.GetBoundaryElement(p_node->GetNextBoundaryElementIndex());
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(0),98);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(1),99);
        p_boundary_element = mesh.GetBoundaryElement(p_node->GetNextBoundaryElementIndex());
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(0),99);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(1),0);
    }
    
    void Test3DMeshCrossReference()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        
        ConformingTetrahedralMesh<3,3> mesh;
        
        mesh.ConstructFromMeshReader(mesh_reader);
        
        
        Node<3> *p_node=mesh.GetNode(34);
        
        TS_ASSERT_EQUALS(p_node->GetNumContainingElements(), 10u);
        Element<3,3> *p_element;
        
        p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),22);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),34);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2),33);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(3),10);
        
        p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),22);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),35);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2),33);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(3),34);
        
        for (int i=0; i<9; i++)
        {
            p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        }
        
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),22);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),34);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2),33);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(3),10);
        
        //Now look at a boundary node
        TS_ASSERT_EQUALS(p_node->GetNumBoundaryElements(), 4u);
        const BoundaryElement<2,3> *p_boundary_element;
        p_boundary_element = mesh.GetBoundaryElement(p_node->GetNextBoundaryElementIndex());
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(0),6);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(1),34);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(2),24);
        p_boundary_element = mesh.GetBoundaryElement(p_node->GetNextBoundaryElementIndex());
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(0),6);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(1),30);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(2),34);
        p_boundary_element = mesh.GetBoundaryElement(p_node->GetNextBoundaryElementIndex());
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(0),24);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(1),34);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(2),10);
        p_boundary_element = mesh.GetBoundaryElement(p_node->GetNextBoundaryElementIndex());
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(0),34);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(1),30);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(2),10);
        
        
        
    }
    
    void Test1DSetPoint()
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        
        ConformingTetrahedralMesh<1,1> mesh;
        
        mesh.ConstructFromMeshReader(mesh_reader);
        
        
        const int node_index=3;
        Node<1> *p_node=mesh.GetNode(node_index);
        
        Point<1> point=p_node->GetPoint();
        TS_ASSERT_DELTA(point[0],0.3,1e-6);
        
        Element<1,1> *p_element;
        p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.1, 1e-6);
        p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.1, 1e-6);
        
        
        // Move node 3 from 0.3 (between node 2 at 0.2 and node 4 at 0.4
        point.SetCoordinate(0,0.25);
        mesh.SetNode(node_index, point);
        p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.05, 1e-6);
        p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.15, 1e-6);
        
        // Move node 3 from 0.3 (between node 2 at 0.2 and node 4 at 0.4
        point.SetCoordinate(0,0.201);
        mesh.SetNode(node_index, point);
        p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.001, 1e-6);
        p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.199, 1e-6);
        
        
        // Move node 3 so that one element is empty
        point.SetCoordinate(0,0.200);
        TS_ASSERT_THROWS_ANYTHING(mesh.SetNode(node_index, point));
        
        // Move node 3 so that one element is negative
        point.SetCoordinate(0,0.15);
        TS_ASSERT_THROWS_ANYTHING(mesh.SetNode(node_index, point));
        
        //Move node 3 back (and recover)
        point.SetCoordinate(0,0.3);
        mesh.SetNode(node_index, point);
        p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.1, 1e-6);
        p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.1, 1e-6);
        
    }
    
    
    void Test2DSetPoint()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        
        ConformingTetrahedralMesh<2,2> mesh;
        
        mesh.ConstructFromMeshReader(mesh_reader);
        
        
        const int node_index=234;
        const int boundary_node_index=99;
        Node<2> *p_node=mesh.GetNode(node_index);
        //Just focus on one element
        Element<2,2> *p_element;
        p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        
        Point<2> point=p_node->GetPoint();
        TS_ASSERT_DELTA(point[0], 0.063497248392600097, 1e-6);
        TS_ASSERT_DELTA(point[1], -0.45483180039309123, 1e-6);
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.00907521, 1e-6);
        
        //Nudge
        point.SetCoordinate(0,0.06);
        mesh.SetNode(node_index, point);
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.00861908, 1e-6);
        
        //Nudge
        point.SetCoordinate(0,0.02);
        mesh.SetNode(node_index, point);
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.00340215, 1e-6);
        
        
        //Nudge
        point.SetCoordinate(0,-0.006);
        mesh.SetNode(node_index, point);
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 1.11485e-05, 1e-6);
        
        //Nudge too far
        point.SetCoordinate(0,-0.0065);
        TS_ASSERT_THROWS_ANYTHING(mesh.SetNode(node_index, point));
        
        
        //Put it back
        point.SetCoordinate(0,0.063497248392600097);
        mesh.SetNode(node_index, point);
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.00907521, 1e-6);
        
        //Now try to move a boundary node
        p_node=mesh.GetNode(boundary_node_index);
        Point<2> boundary_point=p_node->GetPoint();
        TS_ASSERT_DELTA(boundary_point[0], 0.99211470130000001, 1e-6);
        TS_ASSERT_DELTA(boundary_point[1], -0.12533323360000001, 1e-6);
        
        const BoundaryElement<1,2> *p_boundary_element;
        p_boundary_element = mesh.GetBoundaryElement(p_node->GetNextBoundaryElementIndex());
        TS_ASSERT_DELTA(p_boundary_element->GetJacobianDeterminant(), 0.0628215, 1e-6);
        boundary_point.SetCoordinate(0, 1.0);
        mesh.SetNode(boundary_node_index, boundary_point);
        TS_ASSERT_DELTA(p_boundary_element->GetJacobianDeterminant(), 0.0645268, 1e-6);
    }
    
    void TestMovingNodesIn3D()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        double reference_volume = mesh.CalculateMeshVolume();
        
        const int interior_node_index=34;
        Node<3> *p_node=mesh.GetNode(interior_node_index);
        //Just focus on one element
        Element<3,3> *p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        BoundaryElement<2,3> *p_boundary_element = mesh.GetBoundaryElement(p_node->GetNextBoundaryElementIndex());
        
        Point<3> point=p_node->GetPoint();
        TS_ASSERT_DELTA(point[0], 1, 1e-6);
        TS_ASSERT_DELTA(point[1], 0.75, 1e-6);
        TS_ASSERT_DELTA(point[2], 0.75, 1e-6);
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.03125, 1e-6);
        TS_ASSERT_DELTA(p_boundary_element->GetJacobianDeterminant(), 0.125, 1e-6);
        
        // Check the mesh volume hasn't changed
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), reference_volume, 1e-6);
        
        //Nudge
        point.SetCoordinate(2,0.9);
        mesh.SetNode(interior_node_index, point);
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.0125, 1e-6);
        TS_ASSERT_DELTA(p_boundary_element->GetJacobianDeterminant(), 0.05, 1e-6);
        
        //Nudge
        point.SetCoordinate(2, 0.999);
        mesh.SetNode(interior_node_index, point);
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.000125, 1e-6);
        TS_ASSERT_DELTA(p_boundary_element->GetJacobianDeterminant(), 0.0005, 1e-6);
        
        //Nudge
        point.SetCoordinate(2,0.99999);
        mesh.SetNode(interior_node_index, point);
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 1.25e-06, 1e-6);
        TS_ASSERT_DELTA(p_boundary_element->GetJacobianDeterminant(), 5.0e-06, 1e-6);
        
        //Nudge too far
        point.SetCoordinate(2,1.0);
        TS_ASSERT_THROWS_ANYTHING(mesh.SetNode(interior_node_index, point));
        
        //Put it back
        point.SetCoordinate(2,0.75);
        mesh.SetNode(interior_node_index, point);
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.03125, 1e-6);
        TS_ASSERT_DELTA(p_boundary_element->GetJacobianDeterminant(), 0.125, 1e-6);
        
        // Find exterior node
        const int exterior_node_index=0;
        p_node=mesh.GetNode(exterior_node_index); // this exterior node is at (0,0,0)
        point=p_node->GetPoint();
        TS_ASSERT_DELTA(point[0], 0, 1e-6);
        TS_ASSERT_DELTA(point[1], 0, 1e-6);
        TS_ASSERT_DELTA(point[2], 0, 1e-6);
        
        // Move exterior node
        point.SetCoordinate(2,-10.0);
        mesh.SetNode(exterior_node_index, point);
        
        // Check mesh volume has changed
        TS_ASSERT(fabs(mesh.CalculateMeshVolume() - reference_volume) > 1e-1);
    }
    
    void Test1DMeshIn2DSetPoint()
    {
        TrianglesMeshReader<1,2> mesh_reader("mesh/test/data/semicircle_outline");
        
        ConformingTetrahedralMesh<1,2> mesh;
        
        mesh.ConstructFromMeshReader(mesh_reader);
        
        
        const int boundary_node_index=50;
        Node<2> *p_node=mesh.GetNode(boundary_node_index);
        //Just focus on one element
        Element<1,2> *p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        BoundaryElement<0,2> *p_boundary_element = mesh.GetBoundaryElement(p_node->GetNextBoundaryElementIndex());
        
        Point<2> point=p_node->GetPoint();
        TS_ASSERT_DELTA(point[0], -1.0, 1e-6);
        TS_ASSERT_DELTA(point[1], 0.0, 1e-6);
        
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.0628215, 1e-6);
        TS_ASSERT_DELTA(p_boundary_element->GetJacobianDeterminant(), 1.0, 1e-6);
        
        //Nudge left
        point.SetCoordinate(0,-1.5);
        mesh.SetNode(boundary_node_index, point);
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.505885, 1e-6);
        TS_ASSERT_DELTA(p_boundary_element->GetJacobianDeterminant(), 1.0, 1e-6);
        
        //Can't nudge right since an element flips chirality
        point.SetCoordinate(0,-0.5);
        TS_ASSERT_THROWS_ANYTHING(mesh.SetNode(boundary_node_index, point));
        
        
        //Put it back
        point.SetCoordinate(0, -1.0);
        mesh.SetNode(boundary_node_index, point);
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.0628215, 1e-6);
        TS_ASSERT_DELTA(p_boundary_element->GetJacobianDeterminant(), 1.0, 1e-6);
        
    }
    
    void Test2DMeshIn3DSetPoint()
    {
        TrianglesMeshReader<2,3> mesh_reader("mesh/test/data/disk_in_3d");
        
        ConformingTetrahedralMesh<2,3> mesh;
        
        mesh.ConstructFromMeshReader(mesh_reader);
        
        
        const int boundary_node_index=99;
        Node<3> *p_node=mesh.GetNode(boundary_node_index);
        //Just focus on one element
        Element<2,3> *p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        BoundaryElement<1,3> *p_boundary_element = mesh.GetBoundaryElement(p_node->GetNextBoundaryElementIndex());
        
        Point<3> point=p_node->GetPoint();
        TS_ASSERT_DELTA(point[0], 0.99211470130000001, 1e-6);
        TS_ASSERT_DELTA(point[1], -0.12533323360000001, 1e-6);
        TS_ASSERT_DELTA(point[2], 0.0, 1e-6);
        
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.0163772, 1e-6);
        TS_ASSERT_DELTA(p_boundary_element->GetJacobianDeterminant(), 0.0628215, 1e-6);
        
        //Nudge above the plane
        point.SetCoordinate(2,1e-2);
        mesh.SetNode(boundary_node_index, point);
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.0164274 , 1e-6);
        TS_ASSERT_DELTA(p_boundary_element->GetJacobianDeterminant(), 0.0636124, 1e-6);
        
        //Nudge it back
        point.SetCoordinate(2,0.0);
        mesh.SetNode(boundary_node_index, point);
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.0163772, 1e-6);
        TS_ASSERT_DELTA(p_boundary_element->GetJacobianDeterminant(), 0.0628215, 1e-6);
        
        //Nudge below the plane
        point.SetCoordinate(2,-1e-2);
        mesh.SetNode(boundary_node_index, point);
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.0164274, 1e-6);
        TS_ASSERT_DELTA(p_boundary_element->GetJacobianDeterminant(), 0.0636124, 1e-6);
        
        
        //Put it back
        point.SetCoordinate(2,0.0);
        mesh.SetNode(boundary_node_index, point);
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.0163772, 1e-6);
        TS_ASSERT_DELTA(p_boundary_element->GetJacobianDeterminant(), 0.0628215, 1e-6);
        
        //Can't nudge to the other side of the circle without changing handedness
        point.SetCoordinate(0,-1.0);
        point.SetCoordinate(2,0.);
        TS_ASSERT_THROWS_ANYTHING(mesh.SetNode(boundary_node_index, point));
    }
    
    void TestDeletingNodes()
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        Node<1> *p_old_rhs_node = mesh.GetNode(10);
        Node<1> *p_old_lhs_node = mesh.GetNode(0);
        
        ConformingTetrahedralMesh<1,1>::BoundaryElementIterator b_elt_iter;
        ConformingTetrahedralMesh<1,1>::BoundaryNodeIterator b_node_iter;
        
        // Delete the right end node
        mesh.DeleteBoundaryNodeAt(10);
        
        TS_ASSERT(p_old_rhs_node->IsDeleted());
        // Number of *all* nodes & elements should be unchanged, even though we've deleted one,
        // since this is the size of the vector, not the number of active nodes/elements.
        // Yes, it's confusing.
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 11);
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), 10);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 10);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 9);
        // Check the boundary lists are correct
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 2);
        b_node_iter = mesh.GetBoundaryNodeIteratorBegin();
        TS_ASSERT_EQUALS((*b_node_iter++)->GetIndex(), 0);
        TS_ASSERT_EQUALS((*b_node_iter++)->GetIndex(), 9);
        // NB: New boundary elements not added
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 1);
        b_elt_iter = mesh.GetBoundaryElementIteratorBegin();
        TS_ASSERT_EQUALS((*b_elt_iter)->GetNumNodes(), 1);
        TS_ASSERT_EQUALS((*b_elt_iter++)->GetNode(0)->GetIndex(), 0);
        
        // Check the new boundary node
        Node<1> *p_new_rhs_node = mesh.GetNode(9);
        TS_ASSERT(p_new_rhs_node->IsBoundaryNode());
        TS_ASSERT_EQUALS(p_new_rhs_node->GetNumContainingElements(), 1u);
        
        
        
        // Only allowed to remove boundary nodes
        TS_ASSERT_THROWS_ANYTHING(mesh.DeleteBoundaryNodeAt(5));
        
        // Delete the left end node
        mesh.DeleteBoundaryNodeAt(0);
        
        TS_ASSERT(p_old_lhs_node->IsDeleted());
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 11);
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), 10);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 9);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 8);
        // Check the boundary lists are correct
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 2);
        b_node_iter = mesh.GetBoundaryNodeIteratorBegin();
        TS_ASSERT_EQUALS((*b_node_iter++)->GetIndex(), 9); // Note that the boundary is now
        TS_ASSERT_EQUALS((*b_node_iter++)->GetIndex(), 1); // 'reversed'
        // NB: New boundary elements not added
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 0);
        
        // Check the new boundary node
        Node<1> *p_new_lhs_node = mesh.GetNode(1);
        TS_ASSERT(p_new_lhs_node->IsBoundaryNode());
        TS_ASSERT_EQUALS(p_new_lhs_node->GetNumContainingElements(), 1u);
        
        // Check the deleted element/node vectors
        
    }
    
    
    void TestAddingAndDeletingNodes() throw (Exception)
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Add a node at position 0.01
        Point<1> new_point(0.01);
        Element<1,1>* p_first_element = mesh.GetElement(0);
        
        TS_ASSERT_THROWS_NOTHING(mesh.RefineElement(p_first_element, new_point));
        TS_ASSERT_EQUALS(p_first_element->GetNode(1)->GetIndex(), 11);
        // Check the new element is index 10, by comparing nodes
        TS_ASSERT_EQUALS(mesh.GetElement(10)->GetNode(0),
                         p_first_element->GetNode(1));
                         
        // Delete the last node
        mesh.DeleteBoundaryNodeAt(10);
        
        // Add a node
        Point<1> new_point2(0.55);
        Element<1,1>* p_sixth_element = mesh.GetElement(5);
        
        mesh.RefineElement(p_sixth_element, new_point2);
        
        TS_ASSERT_EQUALS(p_sixth_element->GetNode(1)->GetIndex(), 10);
        // Check the new element is index 9, by comparing nodes
        TS_ASSERT_EQUALS(mesh.GetElement(9)->GetNode(0),
                         p_sixth_element->GetNode(1));
                         
    }
    
    
    
    void Test1DNodeMerger()
    {
    
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        
        ConformingTetrahedralMesh<1,1> mesh;
        
        mesh.ConstructFromMeshReader(mesh_reader);
        
        double length=mesh.CalculateMeshVolume();
        const int node_index=3;
        const int target_index=4;
        const int not_neighbour_index=5;
        
        //Cannot merge node 3 with node 5 since they are not neighbours
        TS_ASSERT_THROWS_ANYTHING(mesh.SetNode(node_index, not_neighbour_index));
        
        //Merge node 3 with node 4
        mesh.SetNode(node_index, target_index);
        
        Element<1,1> *p_element;
        p_element = mesh.GetElement(2);
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.2, 1e-6);
        p_element = mesh.GetElement(3);
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.0, 1e-6);
        
        TS_ASSERT_DELTA(length, mesh.CalculateMeshVolume(), 1e-6);
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), mesh.GetNumElements() + 1);
    }
    
    void Test2DNodeMerger()
    {
    
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        
        ConformingTetrahedralMesh<2,2> mesh;
        
        mesh.ConstructFromMeshReader(mesh_reader);
        
        double area=mesh.CalculateMeshVolume();
        //Node 432 is supported by a non-convex region
        //Node 206 is the sole reflex vertex in the region
        //Node 172 is not feasible since it is neighbour to the reflex
        const int node_index=432;
        const int target_index=206;
        const int not_neighbour_index=204;
        const int not_feasible_index=172;
        
        //Element 309 is shared by the moving node (432), the reflex node (206)
        //the non-feasible node (172) - it will vanish
        //Element 762 is shared by the moving node (432), some other node (205)
        //the non-feasible node (172) - it will increase in size
        TS_ASSERT_DELTA(mesh.GetElement(309)->GetJacobianDeterminant(),
                        0.00753493, 1e-6);
        TS_ASSERT_DELTA(mesh.GetElement(762)->GetJacobianDeterminant(),
                        0.00825652, 1e-6);
        //Cannot merge since they are not neighbours
        TS_ASSERT_THROWS_ANYTHING(mesh.SetNode(node_index, not_neighbour_index));
        
        //Cannot merge since an element goes negative
        //The element 763 shared by moving node (432), reflex node (206) and the
        //other neighbour to the reflex node goes negative
        TS_ASSERT_THROWS_ANYTHING(mesh.SetNode(node_index, not_feasible_index, false));
        //Added "crossReference=false" to stop elements deregistering
        
        mesh.SetNode(node_index, target_index);
        
        
        TS_ASSERT_DELTA(area, mesh.CalculateMeshVolume(), 1e-6);
        TS_ASSERT_DELTA(mesh.GetElement(309)->GetJacobianDeterminant(),
                        0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetElement(762)->GetJacobianDeterminant(),
                        0.0126728, 1e-6);
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), mesh.GetNumElements() + 2);
    }
    
    void Test3DNodeMerger() throw (Exception)
    {
    
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_1626_elements");
        
        ConformingTetrahedralMesh<3,3> mesh;
        
        mesh.ConstructFromMeshReader(mesh_reader);
        
        double volume=mesh.CalculateMeshVolume();
        const int node_index=22; //In the middle
        const int target_index=310;
        const int not_neighbour_index=204;
        const int not_feasible_index=103;
        
        TS_ASSERT_THROWS_ANYTHING(mesh.SetNode(node_index, not_neighbour_index));
        TS_ASSERT_THROWS_ANYTHING(mesh.SetNode(node_index, not_feasible_index, false));
        //Added "crossReference=false" to stop elements deregistering
        
        TS_ASSERT_THROWS_NOTHING( mesh.SetNode(node_index, target_index));
        TS_ASSERT_DELTA(volume, mesh.CalculateMeshVolume(), 1e-6);
        
        //Ten elements share 22 and 310.  See:
        /*   Element      N1    N2    N3    N4
                 510     310   348    22   294
                 645      22   328   310   216
                 753      22   329   310   120
                1164     295   310    22   175
                1217     294   310   175    22
                1251     310   336   328    22
                1254     120   310    22   295
                1357     310   336    22   193
                1365      22   329   216   310
                1484     310   348   193    22
        */
        
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), mesh.GetNumElements() + 10);
    }
    
    
    void Test2DBoundaryNodeMergerChangeArea()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        double area=mesh.CalculateMeshVolume();
        double perim=mesh.CalculateMeshSurface();
        int num_nodes=mesh.GetNumNodes();
        int num_elements=mesh.GetNumElements();
        unsigned num_boundary_elements=mesh.GetNumBoundaryElements();
        const int node_index=19;
        const int target_index=20;
        const int not_boundary_index=400;
        TS_ASSERT_THROWS_ANYTHING(mesh.SetNode(node_index, not_boundary_index));
        mesh.SetNode(node_index, target_index);
        
        
        TS_ASSERT_DELTA(area - mesh.CalculateMeshVolume(), 1.24e-4, 1e-6);
        TS_ASSERT_DELTA(perim - mesh.CalculateMeshSurface(), 6.20e-5, 1e-7);
        TS_ASSERT_EQUALS(num_nodes-mesh.GetNumNodes(), 1);
        TS_ASSERT_EQUALS(num_elements-mesh.GetNumElements(), 1);
        TS_ASSERT_EQUALS(num_boundary_elements-mesh.GetNumBoundaryElements(), 1u);
    }
    
    void Test2DBoundaryNodeMerger()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_800_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        double area=mesh.CalculateMeshVolume();
        double perim=mesh.CalculateMeshSurface();
        int num_nodes=mesh.GetNumNodes();
        int num_elements=mesh.GetNumElements();
        unsigned num_boundary_elements=mesh.GetNumBoundaryElements();
        const int node_index=9;
        const int target_index=10;
        const int not_boundary_index=31;
        TS_ASSERT_THROWS_ANYTHING(mesh.SetNode(node_index, not_boundary_index));
        mesh.SetNode(node_index, target_index);
        
        
        TS_ASSERT_DELTA(area - mesh.CalculateMeshVolume(), 0.00, 1e-6);
        TS_ASSERT_DELTA(perim - mesh.CalculateMeshSurface(), 0.00, 1e-7);
        TS_ASSERT_EQUALS(num_nodes-mesh.GetNumNodes(), 1);
        TS_ASSERT_EQUALS(num_elements-mesh.GetNumElements(), 1);
        TS_ASSERT_EQUALS(num_boundary_elements-mesh.GetNumBoundaryElements(), 1u);
        
        const int corner_index=20;
        const int corner_target_index=19;
        mesh.SetNode(corner_index, corner_target_index);
        
        TS_ASSERT_DELTA(area - mesh.CalculateMeshVolume(), 1.25e-5, 1e-7);
        TS_ASSERT_DELTA(perim - mesh.CalculateMeshSurface(), 2.92893e-3, 1e-7);
        TS_ASSERT_EQUALS(num_nodes-mesh.GetNumNodes(), 2);
        TS_ASSERT_EQUALS(num_elements-mesh.GetNumElements(), 2);
        TS_ASSERT_EQUALS(num_boundary_elements-mesh.GetNumBoundaryElements(), 2u);
        
    }
    
    void Test3DBoundaryNodeMerger()
    {
    
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_1626_elements");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        double volume=mesh.CalculateMeshVolume();
        double surface=mesh.CalculateMeshSurface();
        int num_nodes=mesh.GetNumNodes();
        int num_elements=mesh.GetNumElements();
        int num_boundary_elements=mesh.GetNumBoundaryElements();
        const int node_index=147;
        const int target_index=9;
        //const int not_boundary_index=400;
        
        mesh.SetNode(node_index, target_index);
        
        
        TS_ASSERT_DELTA(volume, mesh.CalculateMeshVolume(), 1e-7);
        TS_ASSERT_DELTA(surface, mesh.CalculateMeshSurface(), 1e-7);
        TS_ASSERT_EQUALS(num_nodes-mesh.GetNumNodes(), 1);
        TS_ASSERT_EQUALS(num_elements-mesh.GetNumElements(), 3);
        TS_ASSERT_EQUALS(num_boundary_elements-mesh.GetNumBoundaryElements(), 2);
        
        
        //Can't move corner nodes since this forces some zero volume elements which aren't on the shared list...
    }
    
    
    void TestNodePermutation()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_1626_elements");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        double volume=mesh.CalculateMeshVolume();
        double surface=mesh.CalculateMeshSurface();
        
        Node<3>*  p_node0=mesh.GetNode(0);
        Node<3>*  p_node121=mesh.GetNode(121);
        Node<3>*  p_node125=mesh.GetNode(125);
        Node<3>*  p_node273=mesh.GetNode(273);
        
        RandomNumberGenerator rng;
        mesh.PermuteNodes(rng);
        
        TS_ASSERT_EQUALS(mesh.GetNode(  0)->GetIndex(),   0);
        TS_ASSERT_EQUALS(mesh.GetNode(121)->GetIndex(), 121);
        TS_ASSERT_EQUALS(mesh.GetNode(125)->GetIndex(), 125);
        TS_ASSERT_EQUALS(mesh.GetNode(273)->GetIndex(), 273);
        
        TS_ASSERT_EQUALS(p_node0->GetIndex(), 357);
        TS_ASSERT_EQUALS(p_node121->GetIndex(), 35);
        TS_ASSERT_EQUALS(p_node125->GetIndex(), 219);
        TS_ASSERT_EQUALS(p_node273->GetIndex(), 319);
        TS_ASSERT_EQUALS(mesh.GetNode(p_node0->GetIndex()), p_node0);
        TS_ASSERT_EQUALS(mesh.GetNode(p_node121->GetIndex()), p_node121);
        TS_ASSERT_EQUALS(mesh.GetNode(p_node125->GetIndex()), p_node125);
        TS_ASSERT_EQUALS(mesh.GetNode(p_node273->GetIndex()), p_node273);
        
        TS_ASSERT_DELTA(volume, mesh.CalculateMeshVolume(), 1e-7);
        TS_ASSERT_DELTA(surface, mesh.CalculateMeshSurface(), 1e-7);
      
        
    }
    
        
    void TestConstructRectangle()
    {
        ConformingTetrahedralMesh<2,2> mesh;
        unsigned width=39;
        unsigned height=16;
        mesh.ConstructRectangularMesh(width,height);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), width*height, 1e-7);
        TS_ASSERT_DELTA(mesh.CalculateMeshSurface(), 2.0*(width+height), 1e-7);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(),  (int) ((width+1)*(height+1)));
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(),  (int) (2*(width+height)));
        TS_ASSERT_EQUALS(mesh.GetNumElements(), (int) (2*width*height));
        
        TrianglesMeshWriter<2,2> mesh_writer("","RectangleMesh");
        mesh_writer.WriteFilesUsingMesh(mesh);
    }
   void TestConstructRectangleNoStagger()
    {
        ConformingTetrahedralMesh<2,2> mesh;
        unsigned width=39;
        unsigned height=16;
        mesh.ConstructRectangularMesh(width,height,false);
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), width*height, 1e-7);
        TS_ASSERT_DELTA(mesh.CalculateMeshSurface(), 2.0*(width+height), 1e-7);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(),  (int) ((width+1)*(height+1)));
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(),  (int) (2*(width+height)));
        TS_ASSERT_EQUALS(mesh.GetNumElements(), (int) (2*width*height));
        
        TrianglesMeshWriter<2,2> mesh_writer("","RectangleMeshNoStagger");
        mesh_writer.WriteFilesUsingMesh(mesh);
    }
   
    
    void TestCheckVoronoiDisk()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader,1);
        
        
        TS_ASSERT_EQUALS(mesh.CheckVoronoi(),true);
        
    } 
    
    
    
    void TestCheckVoronoiSquare()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader,1);
        
        
        TS_ASSERT_EQUALS(mesh.CheckVoronoi(),true);
        
    }
    
    void TestCheckCircularFan()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/circular_fan");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader,1);
        
        
        TS_ASSERT_EQUALS(mesh.CheckVoronoi(5e-3),true);
        
    }
    
    void TestCheckMovingMesh()
    {
        ConformingTetrahedralMesh<2,2> mesh;
        
        mesh.ConstructRectangularMesh(1,1);
        
        Node<2> *p_node=mesh.GetNode(1);
        
        Point<2> point=p_node->GetPoint();
        
        for(double x = 1.1; x >= 0.9; x-= 0.01)
        {
            point.SetCoordinate(0,x);
            point.SetCoordinate(1,x);
            mesh.SetNode(1, point);
            
            if(x >= 0.91)
            {
                TS_ASSERT_EQUALS(mesh.CheckVoronoi(0.2),true);
            }
            else
            {
                TS_ASSERT_EQUALS(mesh.CheckVoronoi(0.2),false);
            }
        }
    }
    
    void TestSetOwnerships()
    {
    	TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        unsigned lo = 300;
        unsigned hi = 302;
        
        mesh.SetElementOwnerships(lo, hi);
        
        for (unsigned ele_num=0; ele_num< (unsigned) mesh.GetNumElements(); ele_num++)
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

    void TestConstructCuboid()
    {
        ConformingTetrahedralMesh<3,3> mesh;
        unsigned width=7;
        unsigned height=4;
        unsigned depth=50;
        
        mesh.ConstructCuboid(width,height,depth);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(),  (int) ((width+1)*(height+1)*(depth+1)));
      
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), width*height*depth, 1e-7);
//      TS_ASSERT_DELTA(mesh.CalculateMeshSurface(), 2.0*(width*height+height*depth+depth*width), 1e-7);
        //Each unit square on the surface is split into 2
//      TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(),  (int) (4*(width*height+height*depth+depth*width)));
        //Assuming that each cube is split into 6 tetrahedra
        TS_ASSERT_EQUALS(mesh.GetNumElements(), (int) (6*width*height*depth));

        TrianglesMeshWriter<3,3> mesh_writer("","CuboidMesh");
        mesh_writer.WriteFilesUsingMesh(mesh);
    }
    
};
#endif //_TESTCONFORMINGTETRAHEDRALMESH_HPP_
