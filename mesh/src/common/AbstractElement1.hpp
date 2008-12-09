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


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractElement
{
protected:
    unsigned mIndex;
    std::vector<Node<SPACE_DIM>*> mNodes;

    bool mIsDeleted;
    bool mOwnership;
    bool mFlag;
    unsigned mRegion;
    

    virtual void CommonConstructor(const AbstractElement& rElement)
    {
        mIndex = rElement.mIndex;
        mNodes = rElement.mNodes;

        // Copy various flags
        mIsDeleted = rElement.mIsDeleted;
        mOwnership = rElement.mOwnership;
        mFlag = rElement.mFlag;
        mRegion = rElement.mRegion;
    }

public:
    AbstractElement(unsigned index, const std::vector<Node<SPACE_DIM>*>& rNodes)
        : mIndex(index), mNodes(rNodes)
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
     * Copy constructor. This is needed so that copies of an element don't
     * share pointers to the same matrices, which causes problems when copies
     * get destroyed.
     */
    AbstractElement(const AbstractElement& rElement)
    {
       CommonConstructor(rElement);
    }

    /**
     * \todo Why does the default constructor not do anything?
     */
    AbstractElement()
    {}

    /**
     * Element assignment - make this element equal to the other one.
     */
    virtual AbstractElement& operator=(const AbstractElement& rElement)
    {
        CommonConstructor(rElement);
        return *this;
    }

    virtual ~AbstractElement()
    {}



    virtual void RegisterWithNodes()=0;


    double GetNodeLocation(unsigned localIndex, unsigned dimension) const
    {
        assert(dimension < SPACE_DIM);
        assert((unsigned)localIndex < mNodes.size());
        return mNodes[localIndex]->rGetLocation()[dimension];
    }

    /**
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

    /** Get the index of this element
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
