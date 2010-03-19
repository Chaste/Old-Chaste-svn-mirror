/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
#include "VertexElement.hpp"
#include "RandomNumberGenerator.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>::VertexElement(unsigned index,
												     std::vector<VertexElement<ELEMENT_DIM-1,SPACE_DIM>*> faces,
												     std::vector<bool> orientations)
	  :AbstractElement<ELEMENT_DIM, SPACE_DIM>(index),
      mFaces(faces),
      mOrientations(orientations)
{


    assert(mFaces.size() == mOrientations.size());

    // Populate mNodes using mFaces
    //Make a set of nodes with mFaces
    std::set<Node<SPACE_DIM>* > nodes_set;
    for (unsigned face_index=0; face_index<faces.size(); face_index++)
    {
    	for(unsigned node_index=0; node_index<mFaces[face_index]->GetNumNodes(); node_index++)
     	{
     		nodes_set.insert(mFaces[face_index]->GetNode(node_index));
     	}
    }
    // Populate mNodes
    for (typename std::set< Node<SPACE_DIM>* >::iterator node_iter = nodes_set.begin();
 				node_iter != nodes_set.end();
                ++node_iter)
    {
     	this->mNodes.push_back(*node_iter);
    }
    // Register element with nodes
 	RegisterWithNodes();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>::VertexElement(unsigned index,
													 std::vector<Node<SPACE_DIM>*> nodes)
	: AbstractElement<ELEMENT_DIM, SPACE_DIM>(index, nodes)
{
//	#define COVERAGE_IGNORE
//	assert(SPACE_DIM == 2);
//	#undef COVERAGE_IGNORE

	// \todo this would stop 2d meshes in 3d space
	if(SPACE_DIM == ELEMENT_DIM)
	{
		RegisterWithNodes();
	}
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>::~VertexElement()
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexElement<ELEMENT_DIM, SPACE_DIM>::GetNumFaces() const
{
    return mFaces.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::RegisterWithNodes()
{
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        this->mNodes[i]->AddElement(this->mIndex);
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::MarkAsDeleted()
{
    // Mark element as deleted
    this->mIsDeleted = true;

    // Update nodes in the element so they know they are not contained by it
    for (unsigned i=0; i<this->GetNumNodes(); i++)
    {
        this->mNodes[i]->RemoveElement(this->mIndex);
    }
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::ResetIndex(unsigned index)
{
    for (unsigned i=0; i<this->GetNumNodes(); i++)
    {
       this->mNodes[i]->RemoveElement(this->mIndex);
    }
    this->mIndex = index;
    RegisterWithNodes();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::UpdateNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode)
{
    assert(rIndex < this->mNodes.size());

    // Remove it from the node at this location
    this->mNodes[rIndex]->RemoveElement(this->mIndex);

    // Update the node at this location
    this->mNodes[rIndex] = pNode;

    // Add element to this node
    this->mNodes[rIndex]->AddElement(this->mIndex);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::DeleteNode(const unsigned& rIndex)
{
    assert(rIndex < this->mNodes.size());

    // Remove element from the node at this location
    this->mNodes[rIndex]->RemoveElement(this->mIndex);

    // Remove the node at rIndex (removes node from element)
    this->mNodes.erase(this->mNodes.begin() + rIndex);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::AddNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode)
{
    assert(rIndex < this->mNodes.size());

    // Add pNode to rIndex+1 element of mNodes pushing the others up
    this->mNodes.insert(this->mNodes.begin() + rIndex+1,  pNode);

    // Add element to this node
    this->mNodes[rIndex+1]->AddElement(this->mIndex);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexElement<ELEMENT_DIM, SPACE_DIM>::GetNodeLocalIndex(unsigned globalIndex)
{
    unsigned local_index= UINT_MAX;
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        if (this->GetNodeGlobalIndex(i) == globalIndex)
        {
            local_index = i;
        }
    }
    return local_index;
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class VertexElement<1,1>;
template class VertexElement<1,2>;
template class VertexElement<1,3>;
template class VertexElement<2,2>;
template class VertexElement<2,3>;
template class VertexElement<3,3>;
