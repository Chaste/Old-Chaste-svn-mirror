#ifndef _TESTCONFORMINGTETRAHEDRALMESH_HPP_
#define _TESTCONFORMINGTETRAHEDRALMESH_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"

#include "Element.hpp"
#include "Node.hpp"

#include <vector>

class TestConformingTetrahedralMesh : public CxxTest::TestSuite 
{	
	public:
	
    void TestMeshConstructionFromMeshReader(void)
	{
		TrianglesMeshReader<2,2> meshReader("mesh/test/data/disk_984_elements");
		                  
		ConformingTetrahedralMesh<2,2> mesh;
		
		mesh.ConstructFromMeshReader(meshReader,1);
		
		// Check we have the right number of nodes & elements
		TS_ASSERT_EQUALS(mesh.GetNumCornerNodes(), 543);
		TS_ASSERT_EQUALS(mesh.GetNumElements(), 984);
		
		// Check some node co-ordinates
		TS_ASSERT_DELTA(mesh.GetNodeAt(0)->GetPoint()[0],  0.9980267283, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(0)->GetPoint()[1], -0.0627905195, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(1)->GetPoint()[0], 1.0, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(1)->GetPoint()[1], 0.0, 1e-6);
		
		// Check first element has the right nodes
		ConformingTetrahedralMesh<2,2>::MeshIterator it = mesh.GetElementIteratorBegin();
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(0), 309);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(1), 144);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(2), 310);
		TS_ASSERT_EQUALS(it->GetNode(1), mesh.GetNodeAt(144));
        
        //TS_TRACE("here con tetra\n");
	}
	
	void TestSimpleQuadraticMeshConstructionFromMeshReader(void)
	{
		
		TrianglesMeshReader<2,2> meshReader("mesh/test/data/square_2_elements");
		                  
		ConformingTetrahedralMesh<2,2> mesh;

		try
		{
      		mesh.ConstructFromMeshReader(meshReader,2);
		}
		catch(Exception &e)
		{
			std::cout << e.GetMessage() << std::endl;
			TS_ASSERT(0);
		}
		
		// Check we have the right number of nodes & elements
		TS_ASSERT_EQUALS(mesh.GetNumCornerNodes(), 4);
		TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 9);
		TS_ASSERT_EQUALS(mesh.GetNumElements(), 2);

		// Check some node co-ordinates
		TS_ASSERT_DELTA(mesh.GetNodeAt(0)->GetPoint()[0], 0.0, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(0)->GetPoint()[1], 0.0, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(1)->GetPoint()[0], 1.0, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(1)->GetPoint()[1], 0.0, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(2)->GetPoint()[0], 1.0, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(2)->GetPoint()[1], 1.0, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(3)->GetPoint()[0], 0.0, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(3)->GetPoint()[1], 1.0, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(4)->GetPoint()[0], 0.0, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(4)->GetPoint()[1], 0.5, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(5)->GetPoint()[0], 0.5, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(5)->GetPoint()[1], 0.5, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(6)->GetPoint()[0], 0.5, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(6)->GetPoint()[1], 0.0, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(7)->GetPoint()[0], 1.0, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(7)->GetPoint()[1], 0.5, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(8)->GetPoint()[0], 0.5, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(8)->GetPoint()[1], 1.0, 1e-6);
		
		// Check all elements have the right nodes
		ConformingTetrahedralMesh<2,2>::MeshIterator it = mesh.GetElementIteratorBegin();
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(0), 3);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(1), 0);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(2), 1);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(3), 4);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(4), 5);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(5), 6);
//		TS_ASSERT_EQUALS(it->GetNode(1), mesh.GetNodeAt(144));
        it++;
        TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(0), 1);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(1), 2);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(2), 3);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(3), 7);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(4), 5);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(5), 8);
        //TS_TRACE("here con tetra\n");

	}
	
	void TestQuadraticMeshConstructionFromMeshReader(void)
	{
		
		TrianglesMeshReader<2,2> meshReader("mesh/test/data/disk_984_elements");
		                  
		ConformingTetrahedralMesh<2,2> mesh;

		try
		{
      		mesh.ConstructFromMeshReader(meshReader,2);
		}
		catch(Exception &e)
		{
			std::cout << e.GetMessage() << std::endl;
			TS_ASSERT(0);
		}
		
		// Check we have the right number of nodes & elements
		TS_ASSERT_EQUALS(mesh.GetNumCornerNodes(), 543);
		//TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 543);
		TS_ASSERT_EQUALS(mesh.GetNumElements(), 984);

		// Check some node co-ordinates
		TS_ASSERT_DELTA(mesh.GetNodeAt(0)->GetPoint()[0],  0.9980267283, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(0)->GetPoint()[1], -0.0627905195, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(1)->GetPoint()[0], 1.0, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(1)->GetPoint()[1], 0.0, 1e-6);
		
		// Check first element has the right nodes
		ConformingTetrahedralMesh<2,2>::MeshIterator it = mesh.GetElementIteratorBegin();
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(0), 309);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(1), 144);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(2), 310);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(3), 543);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(4), 544);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(5), 545);
		TS_ASSERT_EQUALS(it->GetNode(1), mesh.GetNodeAt(144));
        it++;
        TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(3), 546);
   		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(4), 547);
        //TS_TRACE("here con tetra\n");

	}
	
	void Test3dLinearMeshConstructionFromMeshReader(void)
	{
		
		TrianglesMeshReader<3,3> meshReader("mesh/test/data/cube_136_elements");
		                  
		//const int DIM = pMeshReader->GetDimension();
		ConformingTetrahedralMesh<3,3> mesh;

		try
		{
      		mesh.ConstructFromMeshReader(meshReader,1);
		}
		catch(Exception &e)
		{
			std::cout << e.GetMessage() << std::endl;
		}
		
		// Check we have the right number of nodes & elements
		TS_ASSERT_EQUALS(mesh.GetNumCornerNodes(), 51);
		//TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 543);
		TS_ASSERT_EQUALS(mesh.GetNumElements(), 136);

		// Check some node co-ordinates
		TS_ASSERT_DELTA(mesh.GetNodeAt(0)->GetPoint()[0], 0.0, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(0)->GetPoint()[1], 0.0, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(0)->GetPoint()[2], 0.0, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(19)->GetPoint()[0], 0.75, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(19)->GetPoint()[1], 0.25, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(19)->GetPoint()[2], 0.0, 1e-6);
		
	}
	
	void Test3dQuadraticMeshConstructionFromMeshReader(void)
	{
		
		TrianglesMeshReader<3,3> meshReader("mesh/test/data/cube_136_elements");
		                  
		//const int DIM = pMeshReader->GetDimension();
		ConformingTetrahedralMesh<3,3> mesh;

		try
		{
      		mesh.ConstructFromMeshReader(meshReader,2);
		}
		catch(Exception &e)
		{
			std::cout << e.GetMessage() << std::endl;
		}
		
		// Check we have the right number of nodes & elements
		TS_ASSERT_EQUALS(mesh.GetNumCornerNodes(), 51);
		//TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 543);
		TS_ASSERT_EQUALS(mesh.GetNumElements(), 136);

		// Check some node co-ordinates
		TS_ASSERT_DELTA(mesh.GetNodeAt(0)->GetPoint()[0], 0.0, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(0)->GetPoint()[1], 0.0, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(0)->GetPoint()[2], 0.0, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(19)->GetPoint()[0], 0.75, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(19)->GetPoint()[1], 0.25, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(19)->GetPoint()[2], 0.0, 1e-6);
		
		// Check first element has the right nodes
		ConformingTetrahedralMesh<3,3>::MeshIterator it = mesh.GetElementIteratorBegin();
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(0), 17);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(1), 10);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(2), 16);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(3), 18);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(4), 51);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(5), 52);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(6), 53);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(7), 54);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(8), 55);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(9), 56);
		TS_ASSERT_EQUALS(it->GetNode(5), mesh.GetNodeAt(52));
        it++;
        TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(5), 58);
   		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(6), 59);
        //TS_TRACE("here con tetra\n");

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
        for(int i=0; i < 11; i++)
        {
            TS_ASSERT_DELTA(mesh.GetNodeAt(i)->GetPoint()[0], 1.5*(i/10.0) , 0.001);
        }
    }
    
    void Test1DClosedMeshIn2DSpace()
    {        
        TrianglesMeshReader<1,2> mesh_reader("mesh/test/data/circle_outline");;
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
        TrianglesMeshReader<1,2> mesh_reader("mesh/test/data/semicircle_outline");;
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
        TrianglesMeshReader<1,1> meshReader("mesh/test/data/1D_0_to_1_10_elements");
                          
        ConformingTetrahedralMesh<1,1> mesh;

        mesh.ConstructFromMeshReader(meshReader);
        
        Node<1> *p_node=mesh.GetNodeAt(0);
        TS_ASSERT_EQUALS(p_node->GetNumContainingElements(), 1);
        unsigned element_index= p_node->GetNextContainingElementIndex();
        TS_ASSERT_EQUALS(element_index, 0);
                  
        Element<1,1> * p_element = mesh.GetElement(element_index);   
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),0);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),1);
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.1, 1e-7);

        Node<1> *p_node2=mesh.GetNodeAt(1);
        TS_ASSERT_EQUALS(p_node2->GetNumContainingElements(), 2);
               
        p_element = mesh.GetElement(p_node2->GetNextContainingElementIndex());
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),0);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),1);  
        
        p_element =mesh.GetElement(p_node2->GetNextContainingElementIndex());
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),1);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),2);
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.1, 1e-7);
        
        // This should wrap back to 1st element
        p_element = mesh.GetElement(p_node2->GetNextContainingElementIndex());
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),0);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),1);    
    }              
    void Test2DMeshCrossReference()
    {
        TrianglesMeshReader<2,2> meshReader("mesh/test/data/disk_984_elements");
                          
        ConformingTetrahedralMesh<2,2> mesh;

        mesh.ConstructFromMeshReader(meshReader);
        
        Node<2> *p_node=mesh.GetNodeAt(234);
        TS_ASSERT_EQUALS(p_node->GetNumContainingElements(), 5);
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
            
    }              
    void Test3DMeshCrossReference()
    {
        TrianglesMeshReader<3,3> meshReader("mesh/test/data/cube_136_elements");
                          
        ConformingTetrahedralMesh<3,3> mesh;

        mesh.ConstructFromMeshReader(meshReader);
        
         
        Node<3> *p_node=mesh.GetNodeAt(34);
 
        TS_ASSERT_EQUALS(p_node->GetNumContainingElements(), 10);
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
        
        for (int i=0; i<9; i++){
            p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        }
         
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(0),22);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(1),34);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(2),33);
        TS_ASSERT_EQUALS(p_element->GetNodeGlobalIndex(3),10);     
    }
    
    void Test1DSetPoint()
    {
        TrianglesMeshReader<1,1> meshReader("mesh/test/data/1D_0_to_1_10_elements");
                          
        ConformingTetrahedralMesh<1,1> mesh;

        mesh.ConstructFromMeshReader(meshReader);
        
        
        Node<1> *p_node=mesh.GetNodeAt(1);
 
        Point<1> point=p_node->GetPoint();
        TS_ASSERT_DELTA(point[0],0.1,1e-7);
 
        Element<1,1> *p_element;
        p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.1, 1e-7);
        p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.1, 1e-7);
        
        
        point.SetCoordinate(0,0.05);
        mesh.SetNode(1, point);
        p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        //\todo TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.05, 1e-7);
        p_element = mesh.GetElement(p_node->GetNextContainingElementIndex());
        //\todo TS_ASSERT_DELTA(p_element->GetJacobianDeterminant(), 0.15, 1e-7);
         
        
    }
                  
};

#endif //_TESTCONFORMINGTETRAHEDRALMESH_HPP_
