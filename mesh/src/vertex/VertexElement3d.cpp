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


#include "VertexElement3d.hpp"


///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

VertexElement3d::VertexElement3d(unsigned index,
                                 std::vector<VertexElement<2,3>*> faces,
                                 std::vector<bool> orientations)
    :  AbstractElement<3, 3>(index),
	   mFaces(faces),
       mOrientations(orientations)
{

    assert(mFaces.size() == mOrientations.size());
    // Populate mNodes using mFaces
    //Make a set of nodes with mFaces
     std::set<Node<3>* > nodes_set;
     for (unsigned face_index=0; face_index<faces.size(); face_index++)
     {
     	for(unsigned node_index=0; node_index<mFaces[face_index]->GetNumNodes(); node_index++)
     	{
     		nodes_set.insert(mFaces[face_index]->GetNode(node_index));
     	}
     }
     // Populate mNodes
     for (std::set<Node<3>* >::iterator node_iter = nodes_set.begin();
 				node_iter != nodes_set.end();
                 ++node_iter)
     {
     	this->mNodes.push_back(*node_iter);
     }
     // Register element with nodes
 	RegisterWithNodes();
}

VertexElement3d::VertexElement3d()
{
}

VertexElement3d::~VertexElement3d()
{
}

unsigned VertexElement3d::GetNumFaces() const
{
    return mFaces.size();
}

VertexElement<2,3>* VertexElement3d::GetFace(unsigned index) const
{
    assert(index < mFaces.size());
	return mFaces[index];
}

bool VertexElement3d::FaceIsOrientatedClockwise(unsigned index) const
{
	assert(index < mOrientations.size());
    return mOrientations[index];
}

void VertexElement3d::RegisterWithNodes()
{
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        this->mNodes[i]->AddElement(this->mIndex);
    }
}

void VertexElement3d::MarkAsDeleted()
{
    // Mark element as deleted
    this->mIsDeleted = true;

    // Update nodes in the element so they know they are not contained by it
    for (unsigned i=0; i<this->GetNumNodes(); i++)
    {
        this->mNodes[i]->RemoveElement(this->mIndex);
    }
}

void VertexElement3d::UpdateNode(const unsigned& rIndex, Node<3>* pNode)
{
    assert(rIndex < this->mNodes.size());

    // Remove it from the node at this location
    this->mNodes[rIndex]->RemoveElement(this->mIndex);

    // Update the node at this location
    this->mNodes[rIndex] = pNode;

    // Add element to this node
    this->mNodes[rIndex]->AddElement(this->mIndex);
}
