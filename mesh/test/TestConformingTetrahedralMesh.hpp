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
	
	static const int DIMENSION=3;
	static const int DIM=2;
	
    
    void TestMeshConstructionFromMeshReader(void)
	{
		TrianglesMeshReader meshReader("mesh/test/data/disk_984_elements");
		                  
		//const int DIM = pMeshReader->GetDimension();
		ConformingTetrahedralMesh<DIM,DIM> mesh;
		
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
		ConformingTetrahedralMesh<DIM,DIM>::MeshIterator it = mesh.GetFirstElement();
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(0), 309);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(1), 144);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(2), 310);
		TS_ASSERT_EQUALS(it->GetNode(1), mesh.GetNodeAt(144));
        
        //TS_TRACE("here con tetra\n");
	}
	
	void TestSimpleQuadraticMeshConstructionFromMeshReader(void)
	{
		
		TrianglesMeshReader meshReader("mesh/test/data/square_2_elements");
		                  
		ConformingTetrahedralMesh<DIM,DIM> mesh;

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
		ConformingTetrahedralMesh<DIM,DIM>::MeshIterator it = mesh.GetFirstElement();
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
		
		TrianglesMeshReader meshReader("mesh/test/data/disk_984_elements");
		                  
		ConformingTetrahedralMesh<DIM,DIM> mesh;

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
	
	void Test3dLinearMeshConstructionFromMeshReader(void)
	{
		
		TrianglesMeshReader meshReader("mesh/test/data/cube_136_elements");
		                  
		//const int DIM = pMeshReader->GetDimension();
		ConformingTetrahedralMesh<DIMENSION,DIMENSION> mesh;

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
		
		TrianglesMeshReader meshReader("mesh/test/data/cube_136_elements");
		                  
		//const int DIM = pMeshReader->GetDimension();
		ConformingTetrahedralMesh<DIMENSION,DIMENSION> mesh;

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
	
	void TestMeshWithBoundaryElements(void)
	{
		TrianglesMeshReader mesh_reader("mesh/test/data/disk_522_elements");
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
        TrianglesMeshReader mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        Point<1> updatedPoint(1.5);
        mesh.RescaleMeshFromBoundaryNode(updatedPoint,10);
        for(int i=0; i < 11; i++)
        {
            TS_ASSERT_DELTA(mesh.GetNodeAt(i)->GetPoint()[0], 1.5*(i/10.0) , 0.001);
        }
    }
	
	// This test is mainly here for performance testing, to check that loading a
	// (relatively) large mesh doesn't take too long.
	void TestLoadingLargeMesh(void)
	{
		TrianglesMeshReader meshReader("mesh/test/data/heart");
		ConformingTetrahedralMesh<3,3> mesh;
		mesh.ConstructFromMeshReader(meshReader, 1);
		
		// Check we have the right number of nodes & elements
		TS_ASSERT_EQUALS(mesh.GetNumCornerNodes(), 63885);
		TS_ASSERT_EQUALS(mesh.GetNumElements(), 322267);
		TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 41812);
	}
};

#endif //_TESTCONFORMINGTETRAHEDRALMESH_HPP_
