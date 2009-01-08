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

#include "AbstractTissue.hpp"
#include "VertexMesh.hpp"

#include <list>

#include <climits> // work around boost bug

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>


/**
 * A facade class encapsulating a vertex-based 'tissue'.
 *
 * Contains a group of cells and maintains the associations 
 * between TissueCells and elements in the VertexMesh.
 * 
 */
template<unsigned DIM>
class VertexBasedTissue : public AbstractTissue<DIM>
{
private:

    /** Vertex-based mesh associated with the tissue. */
    VertexMesh<DIM, DIM>& mrMesh;
    
    /** Results file for elements. */
    out_stream mpElementFile;

public:

    /**
     * Create a new tissue facade from a mesh and collection of cells.
     *
     * There must be precisely one TissueCell for each VertexElement in 
     * the mesh.
     *
     * @param rMesh reference to a VertexMesh
     * @param cells reference to a vector of TissueCells
     * @param validate whether to validate the tissue when it is created (defaults to true)
     */
    VertexBasedTissue(VertexMesh<DIM, DIM>& rMesh,
                      const std::vector<TissueCell>& rCells,
                      bool validate=true);

    /**
     * Destructor, which frees any memory allocated by the constructor.
     */
    virtual ~VertexBasedTissue();

    /**
     * Get a particular VertexElement.
     * 
     * @param elementIndex the global index of the VertexElement
     * 
     * @return a pointer to the VertexElement.
     */
    VertexElement<DIM, DIM>* GetElement(unsigned elementIndex);
    
    /**
     * @return the number of VertexElements in the tissue.
     */
    unsigned GetNumElements();
        
    /**
     * Overridden GetNumNodes() method.
     * 
     * @return the number of nodes in the tissue.
     */    
    unsigned GetNumNodes();
    
    /**
     * Overridden GetNode() method.
     * 
     * @param index global index of the specified node
     * 
     * @return a pointer to the node.
     */
    Node<DIM>* GetNode(unsigned index);
    
    /**
     * Overridden MoveCell() method.
     * 
     * Move a cell to a new location.
     * 
     * @param iter  pointer to the cell to move
     * @param rNewLocation  where to move it to
     */
    void MoveCell(typename AbstractTissue<DIM>::Iterator iter, ChastePoint<DIM>& rNewLocation);
    
    /**
     * Overridden AddCell() method.
     * Add a new cell to the tissue.
     * 
     * @param cell  the cell to add
     * @param newLocation  the position in space at which to put it
     * 
     * @return address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way)
     */
    TissueCell* AddCell(TissueCell cell, c_vector<double,DIM> newLocation);

    /**
     * Remove all cells labelled as dead.
     *
     * Note that after calling this method the tissue will be in an inconsistent state until
     * the equivalent of a 'remesh' is performed! So don't try iterating over cells or anything
     * like that.
     *
     * @return number of cells removed
     */
    unsigned RemoveDeadCells();
    
    /**
     * Remove the VertexElements which have been marked as deleted, perform 
     * any cell rearrangements if required, and update the correspondence 
     * with TissueCells.
     */
    void Update();

    /**
     * Check the consistency of internal data structures. 
     * Each VertexElement must have a TissueCell associated with it.
     */
    virtual void Validate();

    /**
     * Overridden CreateOutputFiles() method.
     * 
     * @param rDirectory  pathname of the output directory, relative to where Chaste output is stored
     * @param rCleanOutputDirectory  whether to delete the contents of the output directory prior to output file creation
     * @param outputCellMutationStates  whether to create a cell mutation state results file
     * @param outputCellTypes  whether to create a cell type results file
     * @param outputCellVariables  whether to create a cell-cycle variable results file
     * @param outputCellCyclePhases  whether to create a cell-cycle phase results file
     * @param outputCellAncestors  whether to create a cell ancestor results file
     */
    void CreateOutputFiles(const std::string &rDirectory,
                           bool rCleanOutputDirectory,
                           bool outputCellMutationStates,
                           bool outputCellTypes,
                           bool outputCellVariables,
                           bool outputCellCyclePhases,
                           bool outputCellAncestors);
    /**
     * Overridden CloseOutputFiles() method.
     * 
     * @param outputCellMutationStates  whether a cell mutation state results file is open
     * @param outputCellTypes  whether a cell type results file is open
     * @param outputCellVariables  whether a cell-cycle variable results file is open
     * @param outputCellCyclePhases  whether a cell-cycle phase results file is open
     * @param outputCellAncestors  whether a cell ancestor results file is open
     */
    void CloseOutputFiles(bool outputCellMutationStates,
                          bool outputCellTypes,
                          bool outputCellVariables,
                          bool outputCellCyclePhases,
                          bool outputCellAncestors);
                          
    /**
     * Overridden WriteResultsToFiles() method.
     * 
     * @param outputCellMutationStates  whether to output cell mutation state results
     * @param outputCellTypes  whether to output cell type results
     * @param outputCellVariables  whether to output cell-cycle variable results
     * @param outputCellCyclePhases  whether to output cell-cycle phase results
     * @param outputCellAncestors  whether to output cell ancestor results
     */
    void WriteResultsToFiles(bool outputCellMutationStates,
                             bool outputCellTypes,
                             bool outputCellVariables,
                             bool outputCellCyclePhases,
                             bool outputCellAncestors);
    
    /**
     * Iterator class allows one to iterate over cells in the tissue.
     * Dereferencing the iterator will give you the current cell.
     * There are also methods to get the VertexElement representing  
     * this cell, and the location of that element.
     */
//    class Iterator
//    {
//    public:
//
//        /**
//         * Dereference the iterator giving you a *reference* to the current cell.
//         * Make sure to use a reference for the result to avoid copying cells unnecessarily.
//         */
//        inline TissueCell& operator*();
//
//        /**
//         * Member access from a pointer.
//         */
//        inline TissueCell* operator->();
//
//        /**
//         * Get a pointer to the element in the mesh which represents this cell.
//         */
//        inline VertexElement<DIM, DIM>* GetElement();
//
//        /**
//         * Get the location in space of this cell.
//         */
//        inline const c_vector<double, DIM>& rGetLocation();
//
//        /**
//         * Comparison not-equal-to.
//         */
//        inline bool operator!=(const Iterator& other);
//
//        /**
//         * Prefix increment operator.
//         */
//        inline Iterator& operator++();
//
//        /**
//         * Constructor for a new iterator.
//         */
//        Iterator(VertexBasedTissue<DIM>& rTissue, std::list<TissueCell>::iterator cellIter);
//
//        /**
//         * The iterator must have a virtual destructor.
//         */
//        virtual ~Iterator()
//        {}
//
//    private:
//
//        /**
//         * Private helper function saying whether we're at the end of the cells.
//         */
//        inline bool IsAtEnd();
//
//        /** The tissue member. */
//        VertexBasedTissue& mrTissue;
//
//        /** Cell iterator member. */
//        std::list<TissueCell>::iterator mCellIter;
//
//        /** Location index member. */
//        unsigned mLocationIndex;
//    };
//
//    /**
//     * @return iterator pointing to the first cell in the tissue
//     */
//    Iterator Begin();
//
//    /**
//     * @return iterator pointing to one past the last cell in the tissue
//     */
//    Iterator End();

};


//////////////////////////////////////////////////////////////////////////////
//                             Iterator class                               //
//////////////////////////////////////////////////////////////////////////////
//
//template<unsigned DIM>
//TissueCell& VertexBasedTissue<DIM>::Iterator::operator*()
//{
//    assert(!IsAtEnd());
//    return *mCellIter;
//}
//
//template<unsigned DIM>
//TissueCell* VertexBasedTissue<DIM>::Iterator::operator->()
//{
//    assert(!IsAtEnd());
//    return &(*mCellIter);
//}
//
//template<unsigned DIM>
//VertexElement<DIM, DIM>* VertexBasedTissue<DIM>::Iterator::GetElement()
//{
//    assert(!IsAtEnd());
//    return mrTissue.GetElement(mLocationIndex);
//}
//
//template<unsigned DIM>
//bool VertexBasedTissue<DIM>::Iterator::operator!=(const VertexBasedTissue<DIM>::Iterator& other)
//{
//    return mCellIter != other.mCellIter;
//}
//
//template<unsigned DIM>
//typename VertexBasedTissue<DIM>::Iterator& VertexBasedTissue<DIM>::Iterator::operator++()
//{
//    ++mCellIter;
//    if (!IsAtEnd())
//    {
//        mLocationIndex = mCellIter->GetLocationIndex();
//    }
//
//    return (*this);
//}
//
//template<unsigned DIM>
//bool VertexBasedTissue<DIM>::Iterator::IsAtEnd()
//{
//    return mCellIter == mrTissue.rGetCells().end();
//}
//
//template<unsigned DIM>
//VertexBasedTissue<DIM>::Iterator::Iterator(VertexBasedTissue& rTissue, std::list<TissueCell>::iterator cellIter)
//    : mrTissue(rTissue),
//      mCellIter(cellIter)
//{
//    // Make sure the tissue isn't empty
//    assert(mrTissue.rGetCells().size() > 0);
//    if (!IsAtEnd())
//    {
//        mLocationIndex = cellIter->GetLocationIndex();
//    }
//}
//
//template<unsigned DIM>
//typename VertexBasedTissue<DIM>::Iterator VertexBasedTissue<DIM>::Begin()
//{
//    return Iterator(*this, this->mCells.begin());
//}
//
//template<unsigned DIM>
//typename VertexBasedTissue<DIM>::Iterator VertexBasedTissue<DIM>::End()
//{
//    return Iterator(*this, this->mCells.end());
//}

#include "TemplatedExport.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VertexBasedTissue)

#endif /*VERTEXBASEDTISSUE_HPP_*/

