#ifndef _MESH_CPP_
#define _MESH_CPP_

#include "Mesh.hpp"
#include "Exception.hpp"

#include <vector>


template<int ELEMENT_DIM, int SPACE_DIM>
Mesh<ELEMENT_DIM, SPACE_DIM>::Mesh()
{
}

template<int ELEMENT_DIM, int SPACE_DIM>
void Mesh<ELEMENT_DIM, SPACE_DIM>::ConstructFromMeshReader(AbstractMeshReader &rMeshReader)
{
	// Check dimension matches the data
	if (SPACE_DIM != rMeshReader.GetDimension())
	{
		throw Exception("Mesh and MeshReader dimensions do not agree.");
	}
	
	// Add nodes
	std::vector<double> coords;
	for (int i=0; i<rMeshReader.GetNumNodes(); i++)
	{
		coords = rMeshReader.GetNextNode();
		mNodes.push_back(Node<SPACE_DIM>(i, Point<SPACE_DIM>(coords), false));
	}
	
	// Add elements
	std::vector<int> node_indices;
	for (int i=0; i<rMeshReader.GetNumElements(); i++)
	{
		node_indices = rMeshReader.GetNextElement();
		std::vector<Node<SPACE_DIM>*> nodes;
		for (int j=0; j<node_indices.size(); j++)
		{
			nodes.push_back(&mNodes[node_indices[j]]);
		}
		mElements.push_back(Element<ELEMENT_DIM,SPACE_DIM>(nodes));
	}
}


template<int ELEMENT_DIM, int SPACE_DIM>
const Node<SPACE_DIM>& Mesh<ELEMENT_DIM, SPACE_DIM>::GetNodeAt(long index) const
{
	assert(index < mNodes.size());
    return mNodes[index];
}

template<int ELEMENT_DIM, int SPACE_DIM>
long Mesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodes()
{
    return mNodes.size();
}

template<int ELEMENT_DIM, int SPACE_DIM>
long Mesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements()
{
    return mElements.size();
}

#endif // _MESH_CPP_
