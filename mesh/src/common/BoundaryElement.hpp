/*
Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _BOUNDARYELEMENT_HPP_
#define _BOUNDARYELEMENT_HPP_

#include "AbstractElement.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class BoundaryElement : public AbstractElement<ELEMENT_DIM, SPACE_DIM>
{

public:
    BoundaryElement(unsigned index, std::vector<Node<SPACE_DIM>*> nodes)
        : AbstractElement<ELEMENT_DIM, SPACE_DIM>(index, nodes)
    {
        RegisterWithNodes();
    }
    
    /**
     * Create a new boundary element from a Node
     * The element has ELEMENT_DIM=0 and
     * SPACE_DIM identical to that of the node from which it is constructed
     * 
     */
    BoundaryElement(unsigned index,
                    Node<SPACE_DIM> *node)
    {
        assert (ELEMENT_DIM == 0);
        
        // Store Node pointer
        this->mNodes.push_back(node);
        
        this->mJacobian(0,0) = 1.0;
        this->mInverseJacobian(0,0) = 1.0;
        this->mWeightedDirection(0) = 1.0;
        this->mJacobianDeterminant = 1.0;
        
        this->mFlag = false;
        this->mIsDeleted=false;
        RegisterWithNodes();
    }
    
    void RegisterWithNodes()
    {
        for (unsigned i=0; i<this->mNodes.size(); i++)
        {
            this->mNodes[i]->AddBoundaryElement(this->mIndex);
        }
    }
    
    void ResetIndex(unsigned index)
    {
        for (unsigned i=0; i<this->GetNumNodes(); i++)
        {
            this->mNodes[i]->RemoveBoundaryElement(this->mIndex);
        }
        this->mIndex=index;
        RegisterWithNodes();
    }
    
    void MarkAsDeleted()
    {
        this->mIsDeleted = true;
        this->mJacobianDeterminant = 0.0;
        // Update nodes in this element so they know they are not contained by us
        for (unsigned i=0; i<this->GetNumNodes(); i++)
        {
            this->mNodes[i]->RemoveBoundaryElement(this->mIndex);
        }
    }
    /** Update node at the given index
     *  @param rIndex is an local index to which node to change
     *  @param pNode is a pointer to the replacement node
     */
    void UpdateNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode)
    {
        assert(rIndex < this->mNodes.size());
        
        // Remove it from the node at this location
        this->mNodes[rIndex]->RemoveBoundaryElement(this->mIndex);
        
        // Update the node at this location
        this->mNodes[rIndex] = pNode;
        
        // Add element to this node
        this->mNodes[rIndex]->AddBoundaryElement(this->mIndex);
    }
    
};



#endif //_BOUNDARYELEMENT_HPP_

