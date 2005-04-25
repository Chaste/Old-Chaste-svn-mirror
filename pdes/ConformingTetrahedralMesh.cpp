#ifndef _CONFORMINGTETRAHEDRALMESH_CPP_
#define _CONFORMINGTETRAHEDRALMESH_CPP_

#include "ConformingTetrahedralMesh.hpp"
#include "Exception.hpp"

#include <vector>


template<int ELEMENT_DIM, int SPACE_DIM>
ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConformingTetrahedralMesh()
{
}

template<int ELEMENT_DIM, int SPACE_DIM>
ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConformingTetrahedralMesh(long numElements)
{
    mElements.reserve(numElements);
}


template<int ELEMENT_DIM, int SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructFromMeshReader(AbstractMeshReader &rMeshReader)
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

//template<int ELEMENT_DIM, int SPACE_DIM>
//ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::~ConformingTetrahedralMesh()
//{
//    delete mpElements;   
//}

template<int ELEMENT_DIM, int SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::AddElement(Element<ELEMENT_DIM, SPACE_DIM> newElement)
{
    mElements.push_back(newElement);
}

template<int ELEMENT_DIM, int SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::AddNode(Node<SPACE_DIM> newNode)
{
	newNode.SetIndex(mNodes.size());
    mNodes.push_back(newNode);
}
//template<int ELEMENT_DIM, int SPACE_DIM>
//void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::AddElement(Element<ELEMENT_DIM, SPACE_DIM> newElement,
//                                                                   std::vector<int> boundaryIndices)
//{
//    mpElements->push_back(newElement);
//    for(int i=0; i < boundaryIndices.size(); i++)
//    {
//        mpBoundaryElements->push_back(mpElements[mpElements->size()]);
//    }
//}
//
//template<int ELEMENT_DIM, int SPACE_DIM>
//void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::AddNode(Node<SPACE_DIM> newNode, bool isBoundaryNode)
//{
//    mpNodes.push_back(newNode);
//    
//}

template<int ELEMENT_DIM, int SPACE_DIM>
const Node<SPACE_DIM>& ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNodeAt(long index) const
{
	assert(index < mNodes.size());
    return mNodes[index];
}

template<int ELEMENT_DIM, int SPACE_DIM>
long ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodes()
{
    return mNodes.size();
}

template<int ELEMENT_DIM, int SPACE_DIM>
long ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements()
{
    return mElements.size();
}

#endif // _CONFORMINGTETRAHEDRALMESH_CPP_
