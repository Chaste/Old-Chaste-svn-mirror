#ifndef _TESTCONFORMINGTETRAHEDRALMESH_HPP_
#define _TESTCONFORMINGTETRAHEDRALMESH_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"


#include "Element.hpp"
#include "Node.hpp"

#include <vector>

class TestConformingTetrahedralMesh : public CxxTest::TestSuite 
{	
	public:
	
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
		std::vector<Node<space_dim>*> nodes;
		// create & add the nodes
		for (int i=0; i < element_dim+1; i++)
		{
			nodes.push_back(new Node<space_dim>(i, false, 0.1+i,0.1+i,0.1+i));
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
			TS_ASSERT_EQUALS(mesh1.GetNodeAt(i).GetIndex(), i);
		}
		for (int i=0; i<mesh2.GetNumNodes(); i++)
		{
			TS_ASSERT_EQUALS(mesh2.GetNodeAt(i).GetIndex(), i);
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
	
};

#endif //_TESTCONFORMINGTETRAHEDRALMESH_HPP_
