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
        TS_ASSERT_EQUALS(mesh.GetNumCornerNodes(), 543U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 984U);
        
        // Check some node co-ordinates
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0],  0.9980267283, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[1], -0.0627905195, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(1)->GetPoint()[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(1)->GetPoint()[1], 0.0, 1e-6);
        
        // Check first element has the right nodes
        ConformingTetrahedralMesh<2,2>::ElementIterator it = mesh.GetElementIteratorBegin();
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(0), 309U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(1), 144U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(2), 310U);
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
        TS_ASSERT_EQUALS(mesh.GetNumCornerNodes(), 4U);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 9U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2U);
        
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
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(0), 3U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(1), 0U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(2), 1U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(3), 4U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(4), 5U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(5), 6U);
//		TS_ASSERT_EQUALS((*it)->GetNode(1), mesh.GetNode(144U));
        it++;
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(0), 1U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(1), 2U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(2), 3U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(3), 7U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(4), 5U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(5), 8U);
        
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
        TS_ASSERT_EQUALS(mesh.GetNumCornerNodes(), 543U);
        //TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 543);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 984U);
        
        // Check some node co-ordinates
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0],  0.9980267283, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[1], -0.0627905195, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(1)->GetPoint()[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(1)->GetPoint()[1], 0.0, 1e-6);
        
        // Check first element has the right nodes
        ConformingTetrahedralMesh<2,2>::ElementIterator it = mesh.GetElementIteratorBegin();
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(0), 309U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(1), 144U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(2), 310U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(3), 543U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(4), 544U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(5), 545U);
        TS_ASSERT_EQUALS((*it)->GetNode(1), mesh.GetNode(144));
        it++;
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(3), 546U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(4), 547U);
        
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
        TS_ASSERT_EQUALS(mesh.GetNumCornerNodes(), 51U);
        //TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 543);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 136U);
        
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
        TS_ASSERT_EQUALS(mesh.GetNumCornerNodes(), 51U);
        //TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 543);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 136U);
        
        // Check some node co-ordinates
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[2], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(19)->GetPoint()[0], 0.75, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(19)->GetPoint()[1], 0.25, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(19)->GetPoint()[2], 0.0, 1e-6);
        
        // Check first element has the right nodes
        ConformingTetrahedralMesh<3,3>::ElementIterator it = mesh.GetElementIteratorBegin();
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(0), 17U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(1), 10U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(2), 16U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(3), 18U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(4), 51U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(5), 52U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(6), 53U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(7), 54U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(8), 55U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(9), 56U);
        TS_ASSERT_EQUALS((*it)->GetNode(5), mesh.GetNode(52));
        it++;
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(5), 58U);
        TS_ASSERT_EQUALS((*it)->GetNodeGlobalIndex(6), 59U);
        
    }
    
    void TestMeshWithBoundaryElements(void)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Check for the right number of boundary edges
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 100U);
        
        // Check all boundary elements have nodes on the boundary
        ConformingTetrahedralMesh<2,2>::BoundaryElementIterator it =
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
        
        TS_ASSERT_EQUALS( mesh.GetNumNodes(), 100U);
        TS_ASSERT_EQUALS( mesh.GetNumElements(), 100U);
        //Check that the mesh_reader has the unculled "faces" (which are nodes)
        TS_ASSERT_EQUALS(mesh_reader.GetNumFaces(), mesh.GetNumNodes());
        TS_ASSERT_EQUALS( mesh.GetNumBoundaryElements(), 0U);
        
    }
    
    void Test1DMeshIn2DSpace()
    {
        TrianglesMeshReader<1,2> mesh_reader("mesh/test/data/semicircle_outline");
        ;
        ConformingTetrahedralMesh<1,2> mesh;
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
        ConformingTetrahedralMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TS_ASSERT_EQUALS( mesh.GetNumNodes(), 132U);
        TS_ASSERT_EQUALS( mesh.GetNumElements(), 224U);
        TS_ASSERT_EQUALS( mesh.GetNumBoundaryElements(), 0U);
    }
    
    void Test2DMeshIn3DSpace()
    {
        TrianglesMeshReader<2,3> mesh_reader("mesh/test/data/disk_in_3d");
        ConformingTetrahedralMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TS_ASSERT_EQUALS( mesh.GetNumNodes(), 312U);
        TS_ASSERT_EQUALS( mesh.GetNumElements(), 522U);
        
        
        ///Check that the mesh_reader has the unculled "faces" (which are edges)
        TS_ASSERT_EQUALS( mesh_reader.GetNumFaces(), 833U);
        //These are the 100 edges around the perimeter of the circle
        TS_ASSERT_EQUALS( mesh.GetNumBoundaryElements(), 100U);
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
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),0U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),1U);
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.1, 1e-6);
        
        Node<1> *p_node2=mesh.GetNode(1);
        TS_ASSERT_EQUALS(p_node2->GetNumContainingElements(), 2u);
        TS_ASSERT_EQUALS(p_node2->GetNumBoundaryElements(), 0u);
        
        p_element = mesh.GetElement(p_node2->GetNextContainingElementIndex());
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),0U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),1U);
        
        p_element =mesh.GetElement(p_node2->GetNextContainingElementIndex());
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),1U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),2U);
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.1, 1e-6);
        
        // This should wrap back to 1st element
        p_element = mesh.GetElement(p_node2->GetNextContainingElementIndex());
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),0U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),1U);
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
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),474U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),290U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2),234U);
        
        p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),234U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),461U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2),460U);
        
        p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),290U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),459U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2),234U);
        
        p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),459U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),461U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2),234U);
        
        p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),460U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),474U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2),234U);
        
        //Now look at a boundary node
        p_node=mesh.GetNode(99);
        TS_ASSERT_EQUALS(p_node->GetNumContainingElements(), 3u);
        TS_ASSERT_EQUALS(p_node->GetNumBoundaryElements(), 2u);
        const BoundaryElement<1,2> *p_boundary_element;
        p_boundary_element = mesh.GetBoundaryElement(p_node->GetNextBoundaryElementIndex());
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(0),98U);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(1),99U);
        p_boundary_element = mesh.GetBoundaryElement(p_node->GetNextBoundaryElementIndex());
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(0),99U);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(1),0U);
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
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),22U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),34U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2),33U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(3),10U);
        
        p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),22U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),35U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2),33U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(3),34U);
        
        for (int i=0; i<9; i++)
        {
            p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        }
        
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),22U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),34U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2),33U);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(3),10U);
        
        //Now look at a boundary node
        TS_ASSERT_EQUALS(p_node->GetNumBoundaryElements(), 4u);
        const BoundaryElement<2,3> *p_boundary_element;
        p_boundary_element = mesh.GetBoundaryElement(p_node->GetNextBoundaryElementIndex());
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(0),6U);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(1),34U);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(2),24U);
        p_boundary_element = mesh.GetBoundaryElement(p_node->GetNextBoundaryElementIndex());
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(0),6U);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(1),30U);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(2),34U);
        p_boundary_element = mesh.GetBoundaryElement(p_node->GetNextBoundaryElementIndex());
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(0),24U);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(1),34U);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(2),10U);
        p_boundary_element = mesh.GetBoundaryElement(p_node->GetNextBoundaryElementIndex());
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(0),34U);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(1),30U);
        TS_ASSERT_EQUALS(p_boundary_element->GetNodeGlobalIndex(2),10U);
        
        
        
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
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 11U);
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), 10U);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 10U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 9U);
        // Check the boundary lists are correct
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 2U);
        b_node_iter = mesh.GetBoundaryNodeIteratorBegin();
        TS_ASSERT_EQUALS((*b_node_iter++)->GetIndex(), 0U);
        TS_ASSERT_EQUALS((*b_node_iter++)->GetIndex(), 9U);
        // NB: New boundary elements not added
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 1U);
        b_elt_iter = mesh.GetBoundaryElementIteratorBegin();
        TS_ASSERT_EQUALS((*b_elt_iter)->GetNumNodes(), 1U);
        TS_ASSERT_EQUALS((*b_elt_iter++)->GetNode(0)->GetIndex(), 0U);
        
        // Check the new boundary node
        Node<1> *p_new_rhs_node = mesh.GetNode(9);
        TS_ASSERT(p_new_rhs_node->IsBoundaryNode());
        TS_ASSERT_EQUALS(p_new_rhs_node->GetNumContainingElements(), 1u);
        
        
        
        // Only allowed to remove boundary nodes
        TS_ASSERT_THROWS_ANYTHING(mesh.DeleteBoundaryNodeAt(5));
        
        // Delete the left end node
        mesh.DeleteBoundaryNodeAt(0);
        
        TS_ASSERT(p_old_lhs_node->IsDeleted());
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 11U);
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), 10U);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 9U);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 8U);
        // Check the boundary lists are correct
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 2U);
        b_node_iter = mesh.GetBoundaryNodeIteratorBegin();
        TS_ASSERT_EQUALS((*b_node_iter++)->GetIndex(), 9U); // Note that the boundary is now
        TS_ASSERT_EQUALS((*b_node_iter++)->GetIndex(), 1U); // 'reversed'
        // NB: New boundary elements not added
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 0U);
        
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
        TS_ASSERT_EQUALS(p_first_element->GetNode(1)->GetIndex(), 11U);
        // Check the new element is index 10, by comparing nodes
        TS_ASSERT_EQUALS(mesh.GetElement(10)->GetNode(0),
                         p_first_element->GetNode(1));
                         
        // Delete the last node
        mesh.DeleteBoundaryNodeAt(10);
        
        // Add a node
        Point<1> new_point2(0.55);
        Element<1,1>* p_sixth_element = mesh.GetElement(5);
        
        mesh.RefineElement(p_sixth_element, new_point2);
        
        TS_ASSERT_EQUALS(p_sixth_element->GetNode(1)->GetIndex(), 10U);
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
        TS_ASSERT_THROWS_ANYTHING(mesh.MoveMergeNode(node_index, not_neighbour_index));
        
        //Merge node 3 with node 4
        mesh.MoveMergeNode(node_index, target_index);
        
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
        TS_ASSERT_THROWS_ANYTHING(mesh.MoveMergeNode(node_index, not_neighbour_index));
        
        //Cannot merge since an element goes negative
        //The element 763 shared by moving node (432), reflex node (206) and the
        //other neighbour to the reflex node goes negative
        TS_ASSERT_THROWS_ANYTHING(mesh.MoveMergeNode(node_index, not_feasible_index, false));
        //Added "crossReference=false" to stop elements deregistering
        
        mesh.MoveMergeNode(node_index, target_index);
        
        
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
        
        TS_ASSERT_THROWS_ANYTHING(mesh.MoveMergeNode(node_index, not_neighbour_index));
        TS_ASSERT_THROWS_ANYTHING(mesh.MoveMergeNode(node_index, not_feasible_index, false));
        //Added "crossReference=false" to stop elements deregistering
        
        TS_ASSERT_THROWS_NOTHING( mesh.MoveMergeNode(node_index, target_index));
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
        unsigned num_nodes=mesh.GetNumNodes();
        unsigned num_elements=mesh.GetNumElements();
        unsigned num_boundary_elements=mesh.GetNumBoundaryElements();
        const int node_index=19;
        const int target_index=20;
        const int not_boundary_index=400;
        TS_ASSERT_THROWS_ANYTHING(mesh.MoveMergeNode(node_index, not_boundary_index));
        mesh.MoveMergeNode(node_index, target_index);
        
        
        TS_ASSERT_DELTA(area - mesh.CalculateMeshVolume(), 1.24e-4, 1e-6);
        TS_ASSERT_DELTA(perim - mesh.CalculateMeshSurface(), 6.20e-5, 1e-7);
        TS_ASSERT_EQUALS(num_nodes-mesh.GetNumNodes(), 1U);
        TS_ASSERT_EQUALS(num_elements-mesh.GetNumElements(), 1U);
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
        unsigned num_elements=mesh.GetNumElements();
        unsigned num_boundary_elements=mesh.GetNumBoundaryElements();
        const int node_index=9;
        const int target_index=10;
        const int not_boundary_index=31;
        TS_ASSERT_THROWS_ANYTHING(mesh.MoveMergeNode(node_index, not_boundary_index));
        mesh.MoveMergeNode(node_index, target_index);
        
        
        TS_ASSERT_DELTA(area - mesh.CalculateMeshVolume(), 0.00, 1e-6);
        TS_ASSERT_DELTA(perim - mesh.CalculateMeshSurface(), 0.00, 1e-7);
        TS_ASSERT_EQUALS(num_nodes-mesh.GetNumNodes(), 1U);
        TS_ASSERT_EQUALS(num_elements-mesh.GetNumElements(), 1U);
        TS_ASSERT_EQUALS(num_boundary_elements-mesh.GetNumBoundaryElements(), 1u);
        
        const int corner_index=20;
        const int corner_target_index=19;
        mesh.MoveMergeNode(corner_index, corner_target_index);
        
        TS_ASSERT_DELTA(area - mesh.CalculateMeshVolume(), 1.25e-5, 1e-7);
        TS_ASSERT_DELTA(perim - mesh.CalculateMeshSurface(), 2.92893e-3, 1e-7);
        TS_ASSERT_EQUALS(num_nodes-mesh.GetNumNodes(), 2U);
        TS_ASSERT_EQUALS(num_elements-mesh.GetNumElements(), 2U);
        TS_ASSERT_EQUALS(num_boundary_elements-mesh.GetNumBoundaryElements(), 2u);
        
    }
    
    void Test3DBoundaryNodeMerger()
    {
    
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_1626_elements");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        double volume=mesh.CalculateMeshVolume();
        double surface=mesh.CalculateMeshSurface();
        unsigned num_nodes=mesh.GetNumNodes();
        unsigned num_elements=mesh.GetNumElements();
        unsigned num_boundary_elements=mesh.GetNumBoundaryElements();
        const int node_index=147;
        const int target_index=9;
        //const int not_boundary_index=400;
        
        mesh.MoveMergeNode(node_index, target_index);
        
        
        TS_ASSERT_DELTA(volume, mesh.CalculateMeshVolume(), 1e-7);
        TS_ASSERT_DELTA(surface, mesh.CalculateMeshSurface(), 1e-7);
        TS_ASSERT_EQUALS(num_nodes-mesh.GetNumNodes(), 1U);
        TS_ASSERT_EQUALS(num_elements-mesh.GetNumElements(), 3U);
        TS_ASSERT_EQUALS(num_boundary_elements-mesh.GetNumBoundaryElements(), 2U);
        
        
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
        
        TS_ASSERT_DELTA(volume, mesh.CalculateMeshVolume(), 1e-7);
        TS_ASSERT_DELTA(surface, mesh.CalculateMeshSurface(), 1e-7);
        
        RandomNumberGenerator::Destroy();
    }
    
    
    void TestConstructRectangle()
    {
        ConformingTetrahedralMesh<2,2> mesh;
        unsigned width=39;
        unsigned height=16;

        mesh.ConstructRectangularMesh(width,height);

        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), width*height, 1e-7);
        TS_ASSERT_DELTA(mesh.CalculateMeshSurface(), 2.0*(width+height), 1e-7);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), ((width+1)*(height+1)));
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 2*(width + height)); 
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(),  2*(width+height) );
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2*width*height);
        
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
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), ((width+1)*(height+1)));
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(),  2*(width+height) );
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2*width*height );
        
        TrianglesMeshWriter<2,2> mesh_writer("","RectangleMeshNoStagger");
        mesh_writer.WriteFilesUsingMesh(mesh);
    }
    
    void TestConstruct1x1RectangularMesh(void)
    {
        ConformingTetrahedralMesh<2,2> rect_mesh;
        rect_mesh.ConstructRectangularMesh(1, 1, false);
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
        
        for (double x = 1.1; x >= 0.9; x-= 0.01)
        {
            point.SetCoordinate(0,x);
            point.SetCoordinate(1,x);
            mesh.SetNode(1, point);
            
            if (x >= 0.91)
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
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        for (unsigned i=0; i<mesh.GetNumBoundaryElements(); i++)
        {
            BoundaryElement<2,3> *b_element=mesh.GetBoundaryElement(i);
            c_vector<double, 3> normal=*(b_element->pGetWeightedDirection());
            c_vector<double, 3> centroid=b_element->CalculateCentroid();
            Point<3> out(centroid+normal);
            Point<3> in(centroid-normal);
            TS_ASSERT_THROWS_NOTHING(mesh.GetContainingElementIndex(in));
            TS_ASSERT_THROWS_ANYTHING(mesh.GetContainingElementIndex(out));
        }
    }
    
    
    void TestConstructCuboid()
    {
    
        ConformingTetrahedralMesh<3,3> mesh;
        unsigned width=7;
        unsigned height=4;
        unsigned depth=5;
        
        unsigned num_boundary_nodes =   2*( (width+1)*(height+1) + (width+1)*(depth+1) + (depth+1)*(height+1) )
                                      - 4*(width-1 + height-1 + depth-1)
                                      - 16;
        
        mesh.ConstructCuboid(width,height,depth);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), ((width+1)*(height+1)*(depth+1)));
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), num_boundary_nodes); 
 
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), width*height*depth, 1e-7);
        TS_ASSERT_DELTA(mesh.CalculateMeshSurface(), 2.0*(width*height+height*depth+depth*width), 1e-7);
        //Each unit square on the surface is split into 2
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(),  4*(width*height+height*depth+depth*width) );
        //Assuming that each cube is split into 6 tetrahedra
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 6*width*height*depth );
        
        for (unsigned i=0; i<mesh.GetNumBoundaryElements(); i++)
        {
            BoundaryElement<2,3> *b_element=mesh.GetBoundaryElement(i);
            c_vector<double, 3> normal=*(b_element->pGetWeightedDirection());
            c_vector<double, 3> centroid=b_element->CalculateCentroid();
            Point<3> out(centroid+normal);
            Point<3> in(centroid-normal);
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
        ConformingTetrahedralMesh<3,3> mesh;
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
        
        mesh.PermuteNodes(perm);
        
        TS_ASSERT_EQUALS(mesh.GetNode(0)->GetIndex(), 0U);
        TS_ASSERT_EQUALS(mesh.GetNode(7)->GetIndex(), 7U);
        
        //Was node 0
        TS_ASSERT_EQUALS(mesh.GetNode(1)->rGetLocation()[0], 0.0);
        TS_ASSERT_EQUALS(mesh.GetNode(1)->rGetLocation()[1], 0.0);
        TS_ASSERT_EQUALS(mesh.GetNode(1)->rGetLocation()[2], 0.0);
        
        //Was node 1
        TS_ASSERT_EQUALS(mesh.GetNode(2)->rGetLocation()[0], 0.2);
        TS_ASSERT_EQUALS(mesh.GetNode(2)->rGetLocation()[1], 0.0);
        TS_ASSERT_EQUALS(mesh.GetNode(2)->rGetLocation()[2], 0.0);
        
        //Was node 2
        TS_ASSERT_EQUALS(mesh.GetNode(0)->rGetLocation()[0], 0.2);
        TS_ASSERT_EQUALS(mesh.GetNode(0)->rGetLocation()[1], 0.2);
        TS_ASSERT_EQUALS(mesh.GetNode(0)->rGetLocation()[2], 0.0);
        
    }
    
    void TestPermuteWithMetisBinaries()
    {
    
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TS_ASSERT_DELTA(mesh.GetNode(0)->rGetLocation()[0],  0.9980, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(0)->rGetLocation()[1], -0.0627, 1e-4);
        mesh.PermuteNodesWithMetisBinaries();
        TS_ASSERT_DELTA(mesh.GetNode(0)->rGetLocation()[0], -0.5358, 1e-4);
        TS_ASSERT_DELTA(mesh.GetNode(0)->rGetLocation()[1], -0.8443, 1e-4);
        
        
        TrianglesMeshReader<3,3> mesh_reader2("mesh/test/data/3D_0_to_.5mm_1889_elements_irregular");
        ConformingTetrahedralMesh<3,3> mesh2;
        mesh2.ConstructFromMeshReader(mesh_reader2);
        
        TS_ASSERT_DELTA(mesh2.GetNode(0)->rGetLocation()[0], 0.0000, 1e-4);
        TS_ASSERT_DELTA(mesh2.GetNode(0)->rGetLocation()[1], 0.0000, 1e-4);
        TS_ASSERT_DELTA(mesh2.GetNode(0)->rGetLocation()[2], 0.0000, 1e-4);
        mesh2.PermuteNodesWithMetisBinaries();
        TS_ASSERT_DELTA(mesh2.GetNode(0)->rGetLocation()[0], 0.0125, 1e-4);
        TS_ASSERT_DELTA(mesh2.GetNode(0)->rGetLocation()[1], 0.0312, 1e-4);
        TS_ASSERT_DELTA(mesh2.GetNode(0)->rGetLocation()[2], 0.0500, 1e-4);
        
        TrianglesMeshWriter<3,3> mesh_writer("","3D_0_to_.5mm_1889_elements_irregular_metis");
        mesh_writer.WriteFilesUsingMesh(mesh2);
    }
    
    void TestDeleteNodes()
    {
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(2,3);
        
        TS_ASSERT_EQUALS(mesh.CalculateMeshVolume(), 6.0);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 12U);
        
        
        //Delete from interior
        mesh.DeleteNode(7);
        TS_ASSERT_EQUALS(mesh.CalculateMeshVolume(), 6.0);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 11U);
        
        //Delete from edge
        mesh.DeleteNode(5);
        TS_ASSERT_EQUALS(mesh.CalculateMeshVolume(), 6.0);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 10U);
        
        //Delete from corner
        mesh.DeleteNode(2);
        TS_ASSERT_EQUALS(mesh.CalculateMeshVolume(), 5.0);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 9U);
        
        // deleting a deleted node should throw an exception
        TS_ASSERT_THROWS_ANYTHING(mesh.DeleteNode(2));
    }
    
    void TestClear()
    {
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(2,3);
        
        TS_ASSERT_EQUALS(mesh.CalculateMeshVolume(), 6.0);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 12u);
        
        mesh.Clear();
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(),0u);
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(),0u);
        TS_ASSERT_EQUALS(mesh.GetNumCornerNodes(),0u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(),0u);
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(),0u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(),0u);
        TS_ASSERT_EQUALS(mesh.GetNumAllBoundaryElements(),0u);
    }
    
    void TestUnflagAllElements()
    {
        ConformingTetrahedralMesh<2,2> mesh;
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
        ConformingTetrahedralMesh<2,2> mesh;
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
    
    void TestCalculateBoundaryOfFlaggedRegion3D()
    {
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(4,4,4);
        mesh.Translate(-2,-2,-2);
        
        // flag elements in the positive octant
        ConformingTetrahedralMesh<3,3>::ElementIterator iter = mesh.GetElementIteratorBegin();
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
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        c_vector<double, 3> location1 = mesh.GetNode(0)->rGetLocation();
        c_vector<double, 3> location2 = mesh.GetNode(2)->rGetLocation();

        // test a normal distance calculation
        c_vector<double, 3> vector = mesh.GetVectorFromAtoB(location1, location2);
        TS_ASSERT_DELTA(vector[0], 0.2, 1e-7);
        TS_ASSERT_DELTA(vector[1], 0.2, 1e-7);
        TS_ASSERT_DELTA(vector[2], 0.0, 1e-7)
        TS_ASSERT_DELTA(norm_2(vector), sqrt(0.08), 1e-7);
        
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
    
    void TestMeshGetWidthMethod(void)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        
        ConformingTetrahedralMesh<2,2> mesh;
        
        mesh.ConstructFromMeshReader(mesh_reader,1);
        
        double width = mesh.GetWidth(1u);
        double height = mesh.GetWidth(2u);
        
        TS_ASSERT_DELTA(width, 2, 1e-6);
        TS_ASSERT_DELTA(height, 2, 1e-6);
        
    }
    
    void TestMeshAddNodeAndReMeshMethod(void)
    {
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(1, 1, false);
                
        TS_ASSERT_EQUALS(mesh.GetNumNodes(),4u);
        
        TS_ASSERT_DELTA(mesh.GetNode(0u)->rGetLocation()[0], 0.0 ,1e-7);
        TS_ASSERT_DELTA(mesh.GetNode(0u)->rGetLocation()[1], 1.0 ,1e-7);
        TS_ASSERT_DELTA(mesh.GetNode(1u)->rGetLocation()[0], 1.0 ,1e-7);
        TS_ASSERT_DELTA(mesh.GetNode(1u)->rGetLocation()[1], 1.0 ,1e-7);
        TS_ASSERT_DELTA(mesh.GetNode(2u)->rGetLocation()[0], 0.0 ,1e-7);
        TS_ASSERT_DELTA(mesh.GetNode(2u)->rGetLocation()[1], 0.0 ,1e-7);
        TS_ASSERT_DELTA(mesh.GetNode(3u)->rGetLocation()[0], 1.0 ,1e-7);
        TS_ASSERT_DELTA(mesh.GetNode(3u)->rGetLocation()[1], 0.0 ,1e-7);
        
        // test the add node method
        c_vector<double ,2> point;
        point[0] = 2.0;
        point[1] = 0.0;
        Node<2>* p_node = new Node<2>(4u, point);
        unsigned new_index = mesh.AddNode(p_node);
        
        TS_ASSERT_EQUALS(new_index, 4u);
        TS_ASSERT_DELTA(mesh.GetNode(4u)->rGetLocation()[0], 2.0 ,1e-7);
        TS_ASSERT_DELTA(mesh.GetNode(4u)->rGetLocation()[1], 0.0 ,1e-7);
        
        // test the add node and ReMesh method
        
        point[0] = 2.0;
        point[1] = 1.0;
        Node<2>* p_node2 = new Node<2>(5u, point);
        
        NodeMap map(mesh.GetNumNodes());
        new_index = mesh.AddNodeAndReMesh(p_node2, map);
        
        TS_ASSERT_EQUALS(new_index, 5u);
        TS_ASSERT_DELTA(mesh.GetNode(5u)->rGetLocation()[0], 2.0 ,1e-7);
        TS_ASSERT_DELTA(mesh.GetNode(5u)->rGetLocation()[1], 1.0 ,1e-7);
        
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 4u);
        
        //TrianglesMeshWriter<2,2> mesh_writer("", "AddNodeAndReMesh", false);
        //mesh_writer.WriteFilesUsingMesh(mesh);
    }
        
};
#endif //_TESTCONFORMINGTETRAHEDRALMESH_HPP_
