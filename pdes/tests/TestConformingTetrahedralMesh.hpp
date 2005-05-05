#ifndef _TESTCONFORMINGTETRAHEDRALMESH_HPP_
#define _TESTCONFORMINGTETRAHEDRALMESH_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.hpp"

#include "Element.hpp"
#include "Node.hpp"

#include <vector>

class TestConformingTetrahedralMesh : public CxxTest::TestSuite 
{	
	public:
	
	/**
	 * \todo This test is obsolete now we use MeshReader. Maybe move some of the
	 * checks into a new test?
	 */
	void TestMeshBuilding(void)
	{
		#define element_dim 3
		#define space_dim 3
		typedef ConformingTetrahedralMesh<element_dim, space_dim>::MeshIterator MeshIterator;
		// Create a mesh
		int n_elements = 5; 
		ConformingTetrahedralMesh<element_dim, space_dim> mesh1(n_elements);
		ConformingTetrahedralMesh<element_dim, space_dim> mesh2;

		std::vector<Element<element_dim, space_dim> > elements;
		std::vector<const Node<space_dim>*> nodes;
		// create & add the nodes
		for (int i=0; i < element_dim+1; i++)
		{
            // Bodge to get positive determinant
            int j=i;
            if (i>1) j=5-i;
			nodes.push_back(new Node<space_dim>(i, false, 0.1+j*j,0.1+j,0.1+j*j*j));
			// Note that the mesh may change the index in general, but
			// shouldn't here since we add in index order.
			mesh1.AddNode(*nodes[i]);
			mesh2.AddNode(*nodes[i]);
		}
		// create the elements and add to the meshes
		for (int i=0; i < n_elements; i++)
		{
			Element<element_dim, space_dim> e(nodes);
			elements.push_back(e);
			mesh1.AddElement(e);
			mesh2.AddElement(e);
		}
		// check nodes have the expected indices
		for (int i=0; i<mesh1.GetNumNodes(); i++)
		{
			TS_ASSERT_EQUALS(mesh1.GetNodeAt(i)->GetIndex(), i);
		}
		for (int i=0; i<mesh2.GetNumNodes(); i++)
		{
			TS_ASSERT_EQUALS(mesh2.GetNodeAt(i)->GetIndex(), i);
		}
		// check elements have the right first node
		// note a typedef (see above) is required for MeshIterator to work.
		for (MeshIterator it=mesh1.GetFirstElement(); it != mesh1.GetLastElement(); it++)
		{
			TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(0), 0);
		}
		for (MeshIterator it=mesh2.GetFirstElement(); it != mesh2.GetLastElement(); it++)
		{
			TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(0), 0);
		}
	}
	
	/**
	 * \todo This test is obsolete now we use MeshReader. Maybe move some of the
	 * checks into a new test?
	 */
	void testNodeMemoryAllocation( void )
    {
        #define DIMENSION 3
        ConformingTetrahedralMesh<DIMENSION,DIMENSION> mesh;
        std::vector< Element<element_dim, space_dim> > elements;
        std::vector< Node<space_dim>* > nodes;
        int n_elements = 5;
        int i = 2;
        const Node<DIMENSION> *a_node, *a_nother_node;
        
        nodes.push_back(new Node<space_dim>(0, false, 0.1+i*i,0.1+i,0.1+i*i*i));
        mesh.AddNode(*nodes[0]);
        a_node = mesh.GetNodeAt(0);
        a_nother_node = mesh.GetNodeAt(0);
        
        TS_ASSERT_EQUALS(a_node,a_nother_node);
        
        
    }
    
    void testMeshConstructionFromMeshReader(void)
	{
		TrianglesMeshReader *pMeshReader = new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_984_elements");
		                  
		//const int DIM = pMeshReader->GetDimension();
		#define DIM 2
		ConformingTetrahedralMesh<DIM,DIM> mesh;
		
		mesh.ConstructFromMeshReader(*pMeshReader,1);
		
		// Check we have the right number of nodes & elements
		TS_ASSERT_EQUALS(mesh.GetNumCornerNodes(), 543);
		TS_ASSERT_EQUALS(mesh.GetNumElements(), 984);
		
		// Check some node co-ordinates
		TS_ASSERT_DELTA(mesh.GetNodeAt(0)->GetPoint()[0],  0.9980267283, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(0)->GetPoint()[1], -0.0627905195, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(1)->GetPoint()[0], 1.0, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(1)->GetPoint()[1], 0.0, 1e-6);
		
		// Check first element has the right nodes
		ConformingTetrahedralMesh<DIM,DIM>::MeshIterator it = mesh.GetFirstElement();
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(0), 309);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(1), 144);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(2), 310);
		TS_ASSERT_EQUALS(it->GetNode(1), mesh.GetNodeAt(144));
        
        //TS_TRACE("here con tetra\n");
	}
	
	void testQuadraticMeshConstructionFromMeshReader(void)
	{
		
		TrianglesMeshReader *pMeshReader = new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_984_elements");
		                  
		//const int DIM = pMeshReader->GetDimension();
		//#define DIM 2
		ConformingTetrahedralMesh<DIM,DIM> mesh;

		try
		{
      		mesh.ConstructFromMeshReader(*pMeshReader,2);
		}
		catch(Exception &e)
		{
			std::cout << e.getMessage() << std::endl;
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
		ConformingTetrahedralMesh<DIM,DIM>::MeshIterator it = mesh.GetFirstElement();
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
	
	void test3dQuadraticMeshConstructionFromMeshReader(void)
	{
		
		TrianglesMeshReader *pMeshReader = new TrianglesMeshReader(
		                  "pdes/tests/meshdata/cube_136_elements");
		                  
		//const int DIM = pMeshReader->GetDimension();
		ConformingTetrahedralMesh<DIMENSION,DIMENSION> mesh;

		try
		{
      		mesh.ConstructFromMeshReader(*pMeshReader,2);
		}
		catch(Exception &e)
		{
			std::cout << e.getMessage() << std::endl;
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
		ConformingTetrahedralMesh<DIMENSION,DIMENSION>::MeshIterator it = mesh.GetFirstElement();
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
	
	void testMeshWithBoundaryElements(void)
	{
		TrianglesMeshReader mesh_reader("pdes/tests/meshdata/disk_522_elements");
		ConformingTetrahedralMesh<2,2> mesh;
		mesh.ConstructFromMeshReader(mesh_reader);
		
		// Check for the right number of boundary edges
		TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 100);
		
		// Check all boundary elements have nodes on the boundary
		ConformingTetrahedralMesh<DIM,DIM>::BoundaryElementIterator it = mesh.GetFirstBoundaryElement();
		while (it != mesh.GetLastBoundaryElement())
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
        TrianglesMeshReader mesh_reader("pdes/tests/meshdata/1D_0_to_1_10_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        Point<1> updatedPoint(1.5);
        mesh.RescaleMeshFromBoundaryNode(updatedPoint,10);
        for(int i=0; i < 11; i++)
        {
            TS_ASSERT_DELTA(mesh.GetNodeAt(i)->GetPoint()[0], 1.5*(i/10.0) , 0.001);
        }
    }
	
};

#endif //_TESTCONFORMINGTETRAHEDRALMESH_HPP_
