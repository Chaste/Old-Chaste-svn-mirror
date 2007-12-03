#ifndef _TESTCONFORMINGTETRAHEDRALMESH_HPP_
#define _TESTCONFORMINGTETRAHEDRALMESH_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include "TrianglesMeshWriter.cpp"
#include "RandomNumberGenerator.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PetscTools.hpp"
#include <cmath>
#include <vector>

class TestConformingTetrahedralMesh : public CxxTest::TestSuite
{
private:

    template<unsigned DIM>
    void EdgeIteratorTest(std::string meshFilename) throw(Exception)
    {
        // create a simple mesh
        TrianglesMeshReader<DIM,DIM> mesh_reader(meshFilename);
        ConformingTetrahedralMesh<DIM,DIM> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // delete one of the nodes and hence an element
        // to check that iterator skips deleted elements
        // this causes element 0 to be deleted which is a good choice for coverage of the begin method
        mesh.DeleteBoundaryNodeAt(0);

        // check that we can iterate over the set of edges
        std::set< std::set< unsigned > > edges_visited;
        
        for (typename ConformingTetrahedralMesh<DIM,DIM>::EdgeIterator edge_iterator=mesh.EdgesBegin();
             edge_iterator!=mesh.EdgesEnd();
             ++edge_iterator)
        {
            std::set<unsigned> node_pair;
            node_pair.insert(edge_iterator.GetNodeA()->GetIndex());
            node_pair.insert(edge_iterator.GetNodeB()->GetIndex());
            
            TS_ASSERT_EQUALS(edges_visited.find(node_pair), edges_visited.end());
            edges_visited.insert(node_pair);
        }
        
        // set up expected node pairs
        std::set< std::set<unsigned> > expected_node_pairs;
        for(unsigned i=0; i<mesh.GetNumAllElements(); i++)
        {
            Element<DIM,DIM>* p_element = mesh.GetElement(i);
            if (!p_element->IsDeleted())
            {
                for(unsigned j=0; j<DIM+1; j++)
                {
                    for(unsigned k=0; k<DIM+1; k++)
                    {
                        unsigned node_A = p_element->GetNodeGlobalIndex(j);
                        unsigned node_B = p_element->GetNodeGlobalIndex(k);
                        
                        if(node_A != node_B)
                        {
                            std::set<unsigned> node_pair;
                            node_pair.insert(node_A);
                            node_pair.insert(node_B);
                            
                            expected_node_pairs.insert(node_pair);
                        }
                    }
                }
            }
        }
        
        TS_ASSERT_EQUALS(edges_visited, expected_node_pairs);       
    }
    
public:

    void TestMeshConstructionFromMeshReader(void)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
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
    

    void Test3dLinearMeshConstructionFromMeshReader(void)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        ConformingTetrahedralMesh<3,3> mesh;
        
        try
        {
            mesh.ConstructFromMeshReader(mesh_reader);
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
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        ChastePoint<1> updatedPoint(1.5);
        mesh.RescaleMeshFromBoundaryNode(updatedPoint,10);
        for (int i=0; i < 11; i++)
        {
            TS_ASSERT_DELTA(mesh.GetNode(i)->GetPoint()[0], 1.5*(i/10.0) , 0.001);
        }
    }
    

    void Test1DClosedMeshIn2DSpace()
    {
        TrianglesMeshReader<1,2> mesh_reader("mesh/test/data/circle_outline");
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
        
        // Check that the mesh_reader has the unculled "faces" (which are edges)
        TS_ASSERT_EQUALS( mesh_reader.GetNumFaces(), 833U);
        // These are the 100 edges around the perimeter of the circle
        TS_ASSERT_EQUALS( mesh.GetNumBoundaryElements(), 100U);
    }
    

    void Test1DMeshCrossReference()
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        ConformingTetrahedralMesh<1,1> mesh;
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
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.1, 1e-6);
        
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
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.1, 1e-6);
        
        // There should be no more containing elements
        TS_ASSERT_EQUALS(++elt_iter, p_node2->ContainingElementsEnd());
    }
    
    
    void Test2DMeshCrossReference()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        ConformingTetrahedralMesh<2,2> mesh;
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
        ConformingTetrahedralMesh<3,3> mesh;
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
    
    
    void Test1DSetPoint()
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        const int node_index=3;
        Node<1> *p_node=mesh.GetNode(node_index);
        
        ChastePoint<1> point=p_node->GetPoint();
        TS_ASSERT_DELTA(point[0],0.3,1e-6);
        
        Element<1,1> *p_element;
        Node<1>::ContainingElementIterator elt_iter = p_node->ContainingElementsBegin();
        p_element = mesh.GetElement(*elt_iter);
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.1, 1e-6);
        p_element = mesh.GetElement(*++elt_iter);
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.1, 1e-6);
        
        // Move node 3 from 0.3 (between node 2 at 0.2 and node 4 at 0.4
        point.SetCoordinate(0,0.25);
        mesh.SetNode(node_index, point);
        elt_iter = p_node->ContainingElementsBegin();
        p_element = mesh.GetElement(*elt_iter);
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.05, 1e-6);
        p_element = mesh.GetElement(*++elt_iter);
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.15, 1e-6);
        
        // Move node 3 from 0.3 (between node 2 at 0.2 and node 4 at 0.4
        point.SetCoordinate(0,0.201);
        mesh.SetNode(node_index, point);
        elt_iter = p_node->ContainingElementsBegin();
        p_element = mesh.GetElement(*elt_iter);
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.001, 1e-6);
        p_element = mesh.GetElement(*++elt_iter);
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
        elt_iter = p_node->ContainingElementsBegin();
        p_element = mesh.GetElement(*elt_iter);
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.1, 1e-6);
        p_element = mesh.GetElement(*++elt_iter);
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
        Node<2>::ContainingElementIterator elt_iter = p_node->ContainingElementsBegin();
        p_element = mesh.GetElement(*elt_iter);
        
        ChastePoint<2> point=p_node->GetPoint();
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
        ChastePoint<2> boundary_point=p_node->GetPoint();
        TS_ASSERT_DELTA(boundary_point[0], 0.99211470130000001, 1e-6);
        TS_ASSERT_DELTA(boundary_point[1], -0.12533323360000001, 1e-6);
        
        Node<2>::ContainingBoundaryElementIterator b_elt_iter = p_node->ContainingBoundaryElementsBegin();
        const BoundaryElement<1,2>* p_boundary_element = mesh.GetBoundaryElement(*b_elt_iter);
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
        Element<3,3> *p_element = mesh.GetElement(*p_node->ContainingElementsBegin());
        BoundaryElement<2,3> *p_boundary_element = mesh.GetBoundaryElement(*p_node->ContainingBoundaryElementsBegin());
        
        ChastePoint<3> point=p_node->GetPoint();
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
        Element<1,2> *p_element = mesh.GetElement(*p_node->ContainingElementsBegin());
        BoundaryElement<0,2> *p_boundary_element = mesh.GetBoundaryElement(*p_node->ContainingBoundaryElementsBegin());
        
        ChastePoint<2> point=p_node->GetPoint();
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
        Element<2,3> *p_element = mesh.GetElement(*p_node->ContainingElementsBegin());
        BoundaryElement<1,3> *p_boundary_element = mesh.GetBoundaryElement(*p_node->ContainingBoundaryElementsBegin());
        
        ChastePoint<3> point=p_node->GetPoint();
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
    
    
    void TestDeleteNodePriorToReMesh() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/circular_fan");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // test it can also delete a boundary node
        mesh.DeleteNodePriorToReMesh(0);
        mesh.DeleteNodePriorToReMesh(11);
        
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(),100u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(),98u);

        NodeMap map(mesh.GetNumNodes());
        mesh.ReMesh(map);
        
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 98u);
    }
    
    
    void TestAddingAndDeletingNodes() throw (Exception)
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Add a node at position 0.01
        ChastePoint<1> new_point(0.01);
        Element<1,1>* p_first_element = mesh.GetElement(0);
        
        TS_ASSERT_THROWS_NOTHING(mesh.RefineElement(p_first_element, new_point));
        TS_ASSERT_EQUALS(p_first_element->GetNode(1)->GetIndex(), 11U);
        // Check the new element is index 10, by comparing nodes
        TS_ASSERT_EQUALS(mesh.GetElement(10)->GetNode(0),
                         p_first_element->GetNode(1));
                         
        // Delete the last node
        mesh.DeleteBoundaryNodeAt(10);
        
        // Add a node
        ChastePoint<1> new_point2(0.55);
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
    
    void TestConstructLine()
    {
        ConformingTetrahedralMesh<1,1> mesh;
        unsigned width=39;;

        mesh.ConstructLinearMesh(width);

        TS_ASSERT_DELTA(mesh.CalculateMeshVolume(), width, 1e-7);
        TS_ASSERT_DELTA(mesh.CalculateMeshSurface(), 0u, 1e-7);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), width+1);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 2u); 
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(),  2u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), width);
        
        TrianglesMeshWriter<1,1> mesh_writer("","LineMesh");
        mesh_writer.WriteFilesUsingMesh(mesh);
    }
    
    
    void TestCheckVoronoiDisk()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.CheckVoronoi(),true);
    }
    

    void TestCheckVoronoiSquare()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TS_ASSERT_EQUALS(mesh.CheckVoronoi(),true);
    }
    
    
    void TestCheckCircularFan()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/circular_fan");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        TS_ASSERT_EQUALS(mesh.CheckVoronoi(5e-3),true);
    }
    
    
    void TestCheckMovingMesh()
    {
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(1,1);
        
        Node<2> *p_node=mesh.GetNode(1);
        ChastePoint<2> point=p_node->GetPoint();
        
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
            ChastePoint<3> out(centroid+normal);
            ChastePoint<3> in(centroid-normal);
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
    
    
    void TestDeleteNodes() throw (Exception)
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
        // moving a deleted node should throw an exception
        TS_ASSERT_THROWS_ANYTHING(mesh.MoveMergeNode(2,1));
    }
    
    
    void TestDeleteNodeFails() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/HalfSquareWithExtraNode");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_THROWS_ANYTHING(mesh.DeleteNode(0));
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
    
    void TestFlagElementsNotContainingNodes()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        ConformingTetrahedralMesh<2,2> mesh;
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
    
    
    void TestMeshGetWidthAndWidthExtremesMethod(void)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        
        ConformingTetrahedralMesh<2,2> mesh;
        
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
        new_index = mesh.AddNode(p_node2);
        mesh.ReMesh(map);
        TS_ASSERT_EQUALS(new_index, 5u);
        TS_ASSERT_DELTA(mesh.GetNode(5u)->rGetLocation()[0], 2.0 ,1e-7);
        TS_ASSERT_DELTA(mesh.GetNode(5u)->rGetLocation()[1], 1.0 ,1e-7);
        
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 4u);
    }
    
    
    void TestReindex()
    {
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(10, 10);
        
        unsigned num_old_nodes = mesh.GetNumNodes();
        
        mesh.DeleteNode(50);
        mesh.DeleteNode(0);
        
        NodeMap map(num_old_nodes);
        
        mesh.ReIndex(map);
        
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), (unsigned)(num_old_nodes-2));
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), (unsigned)(num_old_nodes-2));
        
        TS_ASSERT_EQUALS(map.Size(), num_old_nodes);
        TS_ASSERT_EQUALS(map.IsDeleted(50), true);
        TS_ASSERT_EQUALS(map.IsDeleted(0), true);
        
        for(unsigned i=1; i<num_old_nodes; i++)
        {
            if(i<50)
            {
                TS_ASSERT_EQUALS(map.GetNewIndex(i), (unsigned)(i-1));
            }
            if(i>50)
            {
                TS_ASSERT_EQUALS(map.GetNewIndex(i), (unsigned)(i-2));
            }
        }
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
        ConformingTetrahedralMesh<1,1> mesh;
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
        ConformingTetrahedralMesh<2,2> mesh;
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
        
        ChastePoint<3> out_point(0.1, -10., 0.1);
        TS_ASSERT_EQUALS(element3d.IncludesPoint(out_point), false);
        weights=element3d.CalculateInterpolationWeights(out_point);
        TS_ASSERT_LESS_THAN(0.0, weights[0]);
        TS_ASSERT_LESS_THAN(0.0, weights[1]);
        TS_ASSERT_LESS_THAN(weights[2], 0.0);
        TS_ASSERT_LESS_THAN(0.0, weights[3]);
        
        delete nodes3d[0];
        delete nodes3d[1];
        delete nodes3d[2];
        delete nodes3d[3];
    }
   
   
    void TestPointinMesh3D(void)
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_0_to_1mm_6000_elements");
        ConformingTetrahedralMesh<3,3> mesh;
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
    }
    
    
    void TestFloatingPointIn3D()
    {
        //There's some weird failing behaviour in the refined mesh test
        //This test duplicates it
        
        ConformingTetrahedralMesh<3,3> mesh;
        
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
    void TestArchiveConformingTetrahedralMesh()
    {
        OutputFileHandler handler("archive", false);    // do not clean folder
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "conf_mesh.arch";
        
        // Create an ouput archive
        {
            ConformingTetrahedralMesh<2,2>* const p_mesh = new ConformingTetrahedralMesh<2,2>;
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
            ConformingTetrahedralMesh<2,2> mesh;
            p_mesh->ConstructFromMeshReader(mesh_reader);
            TS_ASSERT_EQUALS(p_mesh->GetNumNodes(),4u);
            TS_ASSERT_EQUALS(p_mesh->GetNumElements(),2u);
            
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            
            output_arch << p_mesh;
            delete p_mesh;
        }
        
        {
            ConformingTetrahedralMesh<2,2>* p_mesh2;
            
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
    
    
    void TestEdgeIterator() throw(Exception)
    {
        EdgeIteratorTest<3>("mesh/test/data/cube_2mm_12_elements");
        EdgeIteratorTest<2>("mesh/test/data/square_4_elements");
        EdgeIteratorTest<1>("mesh/test/data/1D_0_to_1_10_elements");
    }
    

    void TestSetupSmasrmMap()
    {
        // create a mesh on [0,2]x[0,2]
        ConformingTetrahedralMesh<2,2> flagged_mesh; 
        flagged_mesh.ConstructRectangularMesh(50,50);
        flagged_mesh.Scale(2.0/50, 2.0/50);
        
        // flag the [0,1]x[0,1] quadrant
        for(unsigned i=0; i<flagged_mesh.GetNumElements(); i++)
        {
            for(unsigned j=0; j<flagged_mesh.GetElement(i)->GetNumNodes(); j++)
            {
                double x = flagged_mesh.GetElement(i)->GetNode(j)->rGetLocation()[0];
                double y = flagged_mesh.GetElement(i)->GetNode(j)->rGetLocation()[1];
                
                if((x<1.0)&&(y<1.0))
                {
                    flagged_mesh.GetElement(i)->Flag();
                }
            }
        }
        
        // create the smasrm map
        flagged_mesh.SetupSmasrmMap();
        std::map<unsigned, unsigned> smasrm_map=flagged_mesh.rGetSmasrmMap();
        
        // map should be surjective with range {0..n-1} where n is number of flagged nodes 
        // domain should be the set of global node indices of the flagged nodes
        
        std::set<unsigned> expected_domain;
        std::set<unsigned> expected_range;
        for (unsigned node_index = 0; node_index < flagged_mesh.GetNumNodes(); node_index++)
        {
            if (flagged_mesh.GetNode(node_index)->IsFlagged(flagged_mesh))
            {
                expected_domain.insert(node_index);
                expected_range.insert(expected_range.size());
            }
        }
        
        for (std::map<unsigned, unsigned>::iterator map_iterator = smasrm_map.begin();
             map_iterator != smasrm_map.end();
             map_iterator++)
        {
            std::set<unsigned>::iterator set_iterator;
            set_iterator = expected_domain.find(map_iterator->first);
            TS_ASSERT_DIFFERS(set_iterator, expected_domain.end());
            expected_domain.erase(set_iterator);
            
            set_iterator = expected_range.find(map_iterator->second);
            TS_ASSERT_DIFFERS(set_iterator, expected_range.end());
            expected_range.erase(set_iterator);
        }
        
        TS_ASSERT_EQUALS(expected_domain.begin(), expected_domain.end());
        TS_ASSERT_EQUALS(expected_range.begin(), expected_range.end());        
    }
    
    void TestGetAngleBetweenNodes()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_DELTA( mesh.GetAngleBetweenNodes(0,1),  0.0,      1e-12);
        TS_ASSERT_DELTA( mesh.GetAngleBetweenNodes(0,2),  M_PI/4,   1e-12);
        TS_ASSERT_DELTA( mesh.GetAngleBetweenNodes(0,3),  M_PI/2,   1e-12);
        TS_ASSERT_DELTA( mesh.GetAngleBetweenNodes(1,0),  M_PI,     1e-12);
        TS_ASSERT_DELTA( mesh.GetAngleBetweenNodes(2,0), -3*M_PI/4, 1e-12);
        
        TS_ASSERT_THROWS_ANYTHING(mesh.GetAngleBetweenNodes(0,0));
    }
};
#endif //_TESTCONFORMINGTETRAHEDRALMESH_HPP_
