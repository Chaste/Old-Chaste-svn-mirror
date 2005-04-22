#ifndef _TESTMESH_HPP_
#define _TESTMESH_HPP_

#include <cxxtest/TestSuite.h>
#include "TrianglesMeshReader.hpp"
#include "Mesh.cpp"

class TestMesh : public CxxTest::TestSuite
{
	public:
	
	void testMeshConstruction(void)
	{
		TrianglesMeshReader *pMeshReader = new TrianglesMeshReader(
		                  "pdes/tests/meshdata/disk_984_elements");
		                  
		//const int DIM = pMeshReader->GetDimension();
		#define DIM 2
		Mesh<DIM,DIM> mesh = Mesh<DIM,DIM>();
		
		mesh.ConstructFromMeshReader(*pMeshReader);
		
		// Check we have the right number of nodes & elements
		TS_ASSERT_EQUALS(mesh.GetNumNodes(), 543);
		TS_ASSERT_EQUALS(mesh.GetNumElements(), 984);
		
		// Check some node co-ordinates
		TS_ASSERT_DELTA(mesh.GetNodeAt(0).GetPoint()[0],  0.9980267283, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(0).GetPoint()[1], -0.0627905195, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(1).GetPoint()[0], 1.0, 1e-6);
		TS_ASSERT_DELTA(mesh.GetNodeAt(1).GetPoint()[1], 0.0, 1e-6);
		
		// Check first element has the right nodes
		Mesh<DIM,DIM>::MeshIterator it = mesh.GetFirstElement();
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(0), 309);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(1), 144);
		TS_ASSERT_EQUALS(it->GetNodeGlobalIndex(2), 310);
		TS_ASSERT_EQUALS(it->GetNode(1), &(mesh.GetNodeAt(144)));
	}
	
};

#endif //_TESTMESHREADERFUNCTIONALITY_HPP_
