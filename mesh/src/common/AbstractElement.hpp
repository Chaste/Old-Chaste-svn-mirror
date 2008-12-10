/*

Copyright (C) University of Oxford, 2008

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

#ifndef ABSTRACTELEMENT_HPP_
#define ABSTRACTELEMENT_HPP_
#include "Node.hpp"
#include "ChastePoint.hpp"
#include "UblasCustomFunctions.hpp"

#include "Exception.hpp"

#include <vector>
#include <cmath>

// When creating an element within a mesh one needs to specify its global index
// If the element is not used within a mesh the following
// constant is used instead.
const unsigned INDEX_IS_NOT_USED=0;

/**
 * Abstract base class for all elements that can occur in meshes.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractElement
{
protected:
    /** The nodes forming this element */
    std::vector<Node<SPACE_DIM>*> mNodes;
    /** The index of this element within the mesh */
    unsigned mIndex;
    /** A region ID */
    unsigned mRegion;

    /**
     * Whether this element has been deleted, and hence it's location in the
     * mesh can be re-used.
     */
    bool mIsDeleted;
    /** Whether the current process owns this element */
    bool mOwnership;
    /** A flag for the use of higher level algorithms */
    bool mFlag;
    
public:
    AbstractElement(unsigned index, const std::vector<Node<SPACE_DIM>*>& rNodes)
        : mNodes(rNodes), mIndex(index)
    {
        // Sanity checking
        assert(ELEMENT_DIM <= SPACE_DIM);
    
        // Initialise flags.
        // This must be done before the Jacobian calculations, or assertions trip.
        mIsDeleted = false;
        mFlag = false;
        mOwnership = true;
        
        mRegion = 0;
    }

    /**
     * \todo Why does the default constructor not do anything?
     */
    AbstractElement()
        : mIndex(INDEX_IS_NOT_USED),
          mRegion(0),
          mIsDeleted(false),
          mOwnership(true),
          mFlag(false)
    {}

    /**
     * Virtual destructor, since this class has virtual methods.
     * Does nothing special.
     */
    virtual ~AbstractElement()
    {}


    /** 
     *  Update node at the given index
     *  @param rIndex is an local index to which node to change
     *  @param pNode is a pointer to the replacement node
     */
    virtual void UpdateNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode)=0;

    /**
     * Replace one of the nodes in this element with another.
     * 
     * @param pOldNode  pointer to the current node
     * @param pNewNode  pointer to the replacement node
     */
    void ReplaceNode(Node<SPACE_DIM>* pOldNode, Node<SPACE_DIM>* pNewNode)
    {
        //assert(pOldNode != pNewNode); /// \todo this will sometimes trip; is it a logic error?
        for (unsigned i=0; i<this->mNodes.size(); i++)
        {
            if (this->mNodes[i]==pOldNode)
            {
                UpdateNode(i,pNewNode);
                return;
            }
        }
        EXCEPTION("You didn't have that node to start with.");
    }

    /**
     * Mark an element as having been removed from the mesh.
     * Also notify nodes in the element that it has been removed.
     */
    virtual void MarkAsDeleted()=0;


    /**
     * Inform all nodes forming this element that they are in this element.
     */
    virtual void RegisterWithNodes()=0;

    /**
     * Get a single component of the location in space of one of the nodes
     * in this element.
     * 
     * @param localIndex  the index of the node to query, in [0,N) where N
     *   is the number of nodes in this element.
     * @param dimension  the spatial dimension to query.
     */
    double GetNodeLocation(unsigned localIndex, unsigned dimension) const
    {
        assert(dimension < SPACE_DIM);
        assert((unsigned)localIndex < mNodes.size());
        return mNodes[localIndex]->rGetLocation()[dimension];
    }

    /**
     * Get the location in space of one of the nodes in this element.
     * 
     * @param localIndex  the index of the node to query, in [0,N) where N
     *   is the number of nodes in this element.
     *
     * \todo this used to return a reference to a c_vector, in which case a
     * weird error arose where it compiled, ran and passed on some machines
     * but failed the tests (bad_size errors) on another machine.
     */
    c_vector<double, SPACE_DIM> GetNodeLocation(unsigned localIndex) const
    {
        assert((unsigned)localIndex < mNodes.size());
        return mNodes[localIndex]->rGetLocation();
    }
    
    unsigned GetNodeGlobalIndex(unsigned localIndex) const
    {
        assert((unsigned)localIndex < mNodes.size());
        return mNodes[localIndex]->GetIndex();
    }

    Node<SPACE_DIM>* GetNode(unsigned localIndex) const
    {
        assert((unsigned)localIndex < mNodes.size());
        return mNodes[localIndex];
    }

    unsigned GetNumNodes() const
    {
        return mNodes.size();
    }


    void AddNode(Node<SPACE_DIM>* node)
    {
        mNodes.push_back(node);
    }

    bool IsDeleted() const
    {
        return mIsDeleted;
    }

    /**
     *  Get the index of this element
     */
    unsigned GetIndex(void) const
    {
        return mIndex;
    }

    void SetIndex(unsigned index)
    {
        mIndex=index;
    }

    bool GetOwnership() const
    {
        return mOwnership;
    }

    void SetOwnership(bool ownership)
    {
        mOwnership=ownership;
    }

    void Flag()
    {
        mFlag = true;
    }

    void Unflag()
    {
        mFlag = false;
    }

    bool IsFlagged() const
    {
        return mFlag;
    }

    void SetRegion(unsigned region)
    {
        mRegion = region;
    }

    unsigned GetRegion()
    {
        return mRegion;
    }


};

#endif /*ABSTRACTELEMENT_HPP_*/
