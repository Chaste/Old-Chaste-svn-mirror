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
#ifndef VERTEXBASEDTISSUE_HPP_
#define VERTEXBASEDTISSUE_HPP_

#include "TissueCell.hpp"
#include "VertexMesh.hpp"

#include <list>

#include <climits> // work around boost bug


/**
 * A facade class encapsulating a vertex-based 'tissue'.
 *
 * Contains a group of cells and maintains the associations 
 * between TissueCells and elements in the VertexMesh.
 * 
 */
template<unsigned DIM>
class VertexBasedTissue
{
protected:

    /** List of cells */
    std::list<TissueCell> mCells;

    /** Map location (node or VertexElement) indices back to cells */
    std::map<unsigned, TissueCell*> mLocationCellMap;

    /** Vertex-based mesh associated with the tissue. */
    VertexMesh<DIM, DIM>& mrMesh;    

public:

    /**
     * Create a new tissue facade from a mesh and collection of cells.
     *
     * There must be precisely one TissueCell for each VertexElement in 
     * the mesh.
     *
     * @param rMesh a VertexMesh
     * @param cells a vector TissueCells
     * @param validate whether to validate the tissue
     */
    VertexBasedTissue(VertexMesh<DIM, DIM>& rMesh,
                      const std::vector<TissueCell>& rCells,
                      bool validate=true);

    /**
     * Destructor.
     */
    ~VertexBasedTissue();

    VertexElement<DIM, DIM>* GetElement(unsigned elementIndex);
    
    unsigned GetNumElements();
    
    std::list<TissueCell>& rGetCells();
    
    unsigned GetNumRealCells();
        
    /**
     * Get mesh associated with tissue.
     */
    VertexMesh<DIM, DIM>& rGetMesh();

    void Update();

    /**
     * Check consistency of our internal data structures. 
     * Each VertexElement must have a TissueCell associated with it.
     */
    virtual void Validate();
    
    /**
     * Iterator class allows one to iterate over cells in the tissue.
     * Dereferencing the iterator will give you the current cell.
     * There are also methods to get the VertexElement representing  
     * this cell, and the location of that element.
     */
    class Iterator
    {
    public:

        /**
         * Dereference the iterator giving you a *reference* to the current cell.
         * Make sure to use a reference for the result to avoid copying cells unnecessarily.
         */
        inline TissueCell& operator*();

        inline TissueCell* operator->();

        /**
         * Get a pointer to the element in the mesh which represents this cell.
         */
        inline VertexElement<DIM, DIM>* GetElement();

        /**
         * Get the location in space of this cell.
         */
        inline const c_vector<double, DIM>& rGetLocation();

        inline bool operator!=(const Iterator& other);

        /**
         * Prefix increment operator.
         */
        inline Iterator& operator++();

        /**
         * Constructor for a new iterator.
         */
        Iterator(VertexBasedTissue<DIM>& rTissue, std::list<TissueCell>::iterator cellIter);

        /**
         * Must have a virtual destructor.
         */
        virtual ~Iterator()
        {}

    private:

        /**
         * Private helper function which tells us if we're pointing at a real cell.
         * Assumes we are within range (i.e. not at End).
         *
         * Real cells are not deleted.
         */
        virtual inline bool IsRealCell();

        /**
         * Private helper function saying whether we're at the end of the cells.
         */
        inline bool IsAtEnd();

        VertexBasedTissue& mrTissue;
        std::list<TissueCell>::iterator mCellIter;
        unsigned mElementIndex;
    };

    /**
     * @return iterator pointing to the first cell in the tissue
     */
    Iterator Begin();

    /**
     * @return iterator pointing to one past the last cell in the tissue
     */
    Iterator End();

};


//////////////////////////////////////////////////////////////////////////////
//                             Iterator class                               //
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
TissueCell& VertexBasedTissue<DIM>::Iterator::operator*()
{
    assert(!IsAtEnd());
    return *mCellIter;
}

template<unsigned DIM>
TissueCell* VertexBasedTissue<DIM>::Iterator::operator->()
{
    assert(!IsAtEnd());
    return &(*mCellIter);
}

template<unsigned DIM>
VertexElement<DIM, DIM>* VertexBasedTissue<DIM>::Iterator::GetElement()
{
    assert(!IsAtEnd());
    return mrTissue.GetElement(mElementIndex);
}

template<unsigned DIM>
bool VertexBasedTissue<DIM>::Iterator::operator!=(const VertexBasedTissue<DIM>::Iterator& other)
{
    return mCellIter != other.mCellIter;
}

template<unsigned DIM>
typename VertexBasedTissue<DIM>::Iterator& VertexBasedTissue<DIM>::Iterator::operator++()
{
    do
    {
        ++mCellIter;
        if (!IsAtEnd())
        {
            mElementIndex = mCellIter->GetLocationIndex();
        }
    }
    while (!IsAtEnd() && !IsRealCell());

    return (*this);
}

template<unsigned DIM>
bool VertexBasedTissue<DIM>::Iterator::IsRealCell()
{
    return true;
}

template<unsigned DIM>
bool VertexBasedTissue<DIM>::Iterator::IsAtEnd()
{
    return mCellIter == mrTissue.rGetCells().end();
}

template<unsigned DIM>
VertexBasedTissue<DIM>::Iterator::Iterator(VertexBasedTissue& rTissue, std::list<TissueCell>::iterator cellIter)
    : mrTissue(rTissue),
      mCellIter(cellIter)
{
    // Make sure the tissue isn't empty
    assert(mrTissue.rGetCells().size() > 0);
    if (!IsAtEnd())
    {
        mElementIndex = cellIter->GetLocationIndex();
    }
    // Make sure we start at a real cell
    if (mCellIter == mrTissue.rGetCells().begin() && !IsRealCell())
    {
        ++(*this);
    }
}

template<unsigned DIM>
typename VertexBasedTissue<DIM>::Iterator VertexBasedTissue<DIM>::Begin()
{
    return Iterator(*this, this->mCells.begin());
}

template<unsigned DIM>
typename VertexBasedTissue<DIM>::Iterator VertexBasedTissue<DIM>::End()
{
    return Iterator(*this, this->mCells.end());
}


#endif /*VERTEXBASEDTISSUE_HPP_*/

