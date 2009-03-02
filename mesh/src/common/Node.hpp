/*

Copyright (C) University of Oxford, 2005-2009

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


#ifndef _NODE_HPP_
#define _NODE_HPP_

//#include "TetrahedralMesh.hpp"
#include "ChastePoint.hpp"
#include <set>

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class TetrahedralMesh;

/**
 * A node in a finite element mesh.
 */
template<unsigned SPACE_DIM>
class Node
{
private:
    unsigned mIndex;
    unsigned mRegion;

    c_vector<double, SPACE_DIM> mLocation;

    bool mIsBoundaryNode;
    bool mIsDeleted;

    // Set of indices of elements containing this node as a vertex
    std::set<unsigned> mElementIndices;
    std::set<unsigned> mBoundaryElementIndices;

    /**
     * Extraction of commonality between the constructors
     */
    void CommonConstructor(unsigned index, bool isBoundaryNode);

public:
    /**
     * There are many ways of creating a node, depending on how you wish to specify it's
     * spatial location.
     */
    Node(unsigned index, ChastePoint<SPACE_DIM> point, bool isBoundaryNode=false);

    Node(unsigned index, std::vector<double> coords, bool isBoundaryNode=false);

    Node(unsigned index, c_vector<double, SPACE_DIM> location, bool isBoundaryNode=false);

    Node(unsigned index, bool isBoundaryNode=false, double v1=0, double v2=0, double v3=0);

    /**
     * Note: setting the point in space is dangerous.
     * Jacobian and JacobianDeterminant of element need to be updated.
     */
    void SetPoint(ChastePoint<SPACE_DIM> point);

    /**
     * This method should only be called during mesh generation.
     */
    void SetIndex(unsigned index);

    void SetAsBoundaryNode(bool value=true);

    ChastePoint<SPACE_DIM> GetPoint() const;

    /**
     * The returned location may not be modified; if you want that functionality use
     * rGetModifiableLocation instead.
     */
    const c_vector<double, SPACE_DIM>& rGetLocation() const;

    /**
     * If you modify the returned location,
     * Jacobian and JacobianDeterminant of elements need to be updated.
     *
     * Don't forget to assign the result of this call to a reference!
     */
    c_vector<double, SPACE_DIM> &rGetModifiableLocation();

    unsigned GetIndex() const;

    bool IsBoundaryNode() const;

    /**
     * Add an element that contains this node.
     *
     * @param index of the element to add.
     */
    void AddElement(unsigned index);

    /**
     * Remove an element that contains this node.
     *
     * @param index of the element to be removed.
     */
    void RemoveElement(unsigned index);

    /**
     * Remove an boundary element that contains this node.
     *
     * @param index of the boundary element to be removed.
     */
    void RemoveBoundaryElement(unsigned index);

    /**
     * Add an boundary element that contains this node.
     *
     * @param index of the element to add.
     */
    void AddBoundaryElement(unsigned index);

    /**
     * Return a set of indices of elements containing this node as a vertex.
     */
    std::set<unsigned> &rGetContainingElementIndices();

    /**
     * Return a set of indices of boundary elements containing this node as a vertex.
     */
    std::set<unsigned> &rGetContainingBoundaryElementIndices();

    unsigned GetNumContainingElements() const;

    unsigned GetNumBoundaryElements() const;

    /**
     * Mark a node as having been removed from the mesh
     */
    void MarkAsDeleted();

    bool IsDeleted() const;

    /**
     * Determine if a node lives within a flagged element.
     */
    template <unsigned ELEMENT_DIM>
    bool IsFlagged(TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh)
    {
        bool in_flagged_element = false;
        for (ContainingElementIterator it = ContainingElementsBegin();
             it != ContainingElementsEnd();
             ++it)
        {
            if (rMesh.GetElement(*it)->IsFlagged())
            {
                in_flagged_element = true;
                break;
            }
        }
        return in_flagged_element;
    }
    
    void SetRegion(unsigned region);

    unsigned GetRegion() const;

    /**
     * An iterator over the indices of elements which contain this node.
     */
    class ContainingElementIterator
    {
    public:
        ContainingElementIterator(std::set<unsigned>::const_iterator indexIterator)
            : mIndexIterator(indexIterator)
        {}

        const unsigned& operator*() const
        {
            return *mIndexIterator;
        }

        bool operator!=(const ContainingElementIterator& other) const
        {
            return mIndexIterator != other.mIndexIterator;
        }
        bool operator==(const ContainingElementIterator& other) const
        {
            return !operator!=(other);
        }

        ContainingElementIterator& operator++()
        {
            ++mIndexIterator;
            return *this;
        }
    private:
        std::set<unsigned>::const_iterator mIndexIterator;
    };

    ContainingElementIterator ContainingElementsBegin() const
    {
        return ContainingElementIterator(mElementIndices.begin());
    }

    ContainingElementIterator ContainingElementsEnd() const
    {
        return ContainingElementIterator(mElementIndices.end());
    }

    /**
     * An iterator over the indices of boundary elements which contain this node.
     */
    class ContainingBoundaryElementIterator
    {
    public:
        ContainingBoundaryElementIterator(std::set<unsigned>::const_iterator indexIterator)
            : mIndexIterator(indexIterator)
        {}

 
        const unsigned& operator*() const
        {
            return *mIndexIterator;
        }

        bool operator!=(const ContainingBoundaryElementIterator& other) const
        {
            return mIndexIterator != other.mIndexIterator;
        }
        bool operator==(const ContainingBoundaryElementIterator& other) const
        {
            return !operator!=(other);
        }

        ContainingBoundaryElementIterator& operator++()
        {
            ++mIndexIterator;
            return *this;
        }
    private:
        std::set<unsigned>::const_iterator mIndexIterator;
    };

    ContainingBoundaryElementIterator ContainingBoundaryElementsBegin() const
    {
        return ContainingBoundaryElementIterator(mBoundaryElementIndices.begin());
    }

    ContainingBoundaryElementIterator ContainingBoundaryElementsEnd() const
    {
        return ContainingBoundaryElementIterator(mBoundaryElementIndices.end());
    }
};


#endif //_NODE_HPP_
