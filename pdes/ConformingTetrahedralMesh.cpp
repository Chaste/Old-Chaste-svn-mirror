#ifndef _CONFORMINGTETRAHEDRALMESH_CPP_
#define _CONFORMINGTETRAHEDRALMESH_CPP_

#include "ConformingTetrahedralMesh.hpp"
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

//template<int ELEMENT_DIM, int SPACE_DIM>
//ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::~ConformingTetrahedralMesh()
//{
//    delete mpElements;   
//}

template<int ELEMENT_DIM, int SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::AddElement(Element<ELEMENT_DIM, SPACE_DIM>& rNewElement)
{
    mElements.push_back(rNewElement);
}

template<int ELEMENT_DIM, int SPACE_DIM>
void ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::AddNode(Node<SPACE_DIM>& rNewNode)
{
	rNewNode.SetIndex(mNodes.size());
    mNodes.push_back(rNewNode);
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
