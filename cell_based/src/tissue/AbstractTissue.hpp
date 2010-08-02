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
#ifndef ABSTRACTTISSUE_HPP_
#define ABSTRACTTISSUE_HPP_

#include "TissueCell.hpp"
#include "OutputFileHandler.hpp"

#include <list>
#include <map>
#include <vector>
#include <boost/shared_ptr.hpp>

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include <boost/serialization/vector.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "Node.hpp"
#include "CellPropertyRegistry.hpp"

// Needed here to avoid serialization errors (on Boost<1.37)
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "ApoptoticCellProperty.hpp"
#include "CellLabel.hpp"

/**
 * An abstract facade class encapsulating a tissue.
 *
 * Contains a group of cells and associated methods.
 *
 */
template<unsigned DIM>
class AbstractTissue
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mCells;
        archive & mLocationCellMap;
        archive & mCellLocationMap;
        archive & mCellProliferativeTypeCount;
        archive & mCellCyclePhaseCount;
        archive & mTissueContainsMesh;
        archive & mpCellPropertyRegistry;
        archive & mOutputCellIdData;
        archive & mOutputCellMutationStates;
        archive & mOutputCellAncestors;
        archive & mOutputCellProliferativeTypes;
        archive & mOutputCellVariables;
        archive & mOutputCellCyclePhases;
        archive & mOutputCellAges;
        archive & mOutputCellVolumes;
    }

protected:

    /** List of cells */
    std::list<TissueCellPtr> mCells;

    /** Map location (node or VertexElement) indices back to cells */
    std::map<unsigned, TissueCellPtr> mLocationCellMap;

    /** Map cells to location (node or VertexElement) indices */
    std::map<TissueCell*, unsigned> mCellLocationMap;

    /** Current cell type counts */
    std::vector<unsigned> mCellProliferativeTypeCount;

    /** Current cell cycle phase counts */
    std::vector<unsigned> mCellCyclePhaseCount;

    /** Results file for node visualization */
    out_stream mpVizNodesFile;

    /** Results file for cell visualization */
    out_stream mpVizCellProliferativeTypesFile;

    /** Results file for cell mutation states */
    out_stream mpCellMutationStatesFile;

    /** Results file for cell ancestors */
    out_stream mpVizCellAncestorsFile;

    /** Results file for cell types */
    out_stream mpCellProliferativeTypesFile;

    /** Results file for cell cycle phases */
    out_stream mpCellCyclePhasesFile;

    /** Results file for cell variables */
    out_stream mpCellVariablesFile;

    /** Results file for cell ages */
    out_stream mpCellAgesFile;

    /** Results file for logged cell data. */
    out_stream mpCellIdFile;

    /** Results file for boundary nodes. */
    out_stream mpVizBoundaryNodesFile;

    /** Whether the tissue contains a mesh. */
    bool mTissueContainsMesh;

    /** Cell property registry. */
    boost::shared_ptr<CellPropertyRegistry> mpCellPropertyRegistry;

    /** Whether to write cell ID data to file. */
    bool mOutputCellIdData;

    /** Whether to count the number of each cell mutation state and output to file. */
    bool mOutputCellMutationStates;

    /** Whether to output the ancestor of each cell to a visualizer file. */
    bool mOutputCellAncestors;

    /** Whether to count the number of each cell type and output to file. */
    bool mOutputCellProliferativeTypes;

    /** Whether to write the cell variables to a file. */
    bool mOutputCellVariables;

    /** Whether to write the cell cycle phases to a file. */
    bool mOutputCellCyclePhases;

    /** Whether to write the cell ages to a file. */
    bool mOutputCellAges;

    /** Whether to write the cell volumes (in 3D) or areas (in 2D) to a file. */
    bool mOutputCellVolumes;

    /**
     * Check consistency of our internal data structures.
     */
    virtual void Validate()=0;
    
    /**
     * Constructor for use by archiving only. Please use the other constructor.
     *
     * Doesn't take in cells, since these are dealt
     * with by the serialize method.
     */
    AbstractTissue()
    {
    }

public:

    /**
     * AbstractTissue Constructor.
     *
     * @note Warning: the passed-in vector of cells will be emptied, even if the constructor
     * throws an exception!
     *
     * @param rCells a vector of cells.  Copies of the cells will be stored in the tissue,
     *     and the passed-in vector cleared.
     * @param locationIndices an optional vector of location indices that correspond to real cells
     */
    AbstractTissue(std::vector<TissueCellPtr>& rCells,
                   const std::vector<unsigned> locationIndices=std::vector<unsigned>());

    /**
     * Base class with virtual methods needs a virtual destructor.
     */
    virtual ~AbstractTissue()
    {}

    /**
     * Initialise each cell's cell cycle model.
     */
    void InitialiseCells();

    /**
     * @return reference to mCells.
     */
    std::list<TissueCellPtr>& rGetCells();

    /**
     * \todo This method returns true if the tissue is a MeshBasedTissue or a
     *       VertexBasedTissue, but is actually used to tell force laws whether
     *       the tissue is a MeshBasedTissue. See #1303.
     *
     * @return whether the tissue contains a mesh.
     */
    bool HasMesh();

    /**
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @return the number of nodes in the tissue.
     */
    virtual unsigned GetNumNodes()=0;

    /**
     * Find where a given cell is in space.
     *
     * @param pCell the cell
     *
     * @return the location of the cell
     */
    virtual c_vector<double, DIM> GetLocationOfCellCentre(TissueCellPtr pCell)=0;

    /**
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param index  global index of the specified node
     *
     * @return a pointer to the node with a given index.
     */
    virtual Node<DIM>* GetNode(unsigned index)=0;

    /**
     * Add a new node to the tissue.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param pNewNode pointer to the new node
     *
     * @return global index of new node in tissue.
     */
    virtual unsigned AddNode(Node<DIM>* pNewNode)=0;

    /**
     * Move the node with a given index to a new point in space.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param nodeIndex the index of the node to be moved
     * @param rNewLocation the new target location of the node
     */
    virtual void SetNode(unsigned nodeIndex, ChastePoint<DIM>& rNewLocation)=0;

    /**
     * Helper method for establishing if a cell is real.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param pCell the cell
     *
     * @return whether a given cell is associated with a deleted
     *         node (cell-centre models) or element (vertex models).
     */
    virtual bool IsCellAssociatedWithADeletedLocation(TissueCellPtr pCell)=0;

    /**
     * Update the location of each node in the tissue given
     * a vector of forces on nodes and a time step over which
     * to integrate the equations of motion.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rNodeForces  forces on nodes
     * @param dt time step
     */
    virtual void UpdateNodeLocations(const std::vector< c_vector<double, DIM> >& rNodeForces, double dt)=0;

    /**
     * Get the damping constant for this node - ie d in drdt = F/d.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param nodeIndex the global index of this node
     *
     * @return the damping constant at the node.
     */
    virtual double GetDampingConstant(unsigned nodeIndex)=0;

    /**
     * Add a new cell to the tissue.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param pNewCell  the cell to add
     * @param rCellDivisionVector  a vector providing information regarding how the cell division should occur
     *     (for cell-centre tissues, this vector is the position of the daughter cell; for vertex tissues it
     *      can be used by any subclass of TissueSimulation to as a means of dictating the axis along which
     *      the parent cell divides)
     * @param pParentCell pointer to a parent cell (if required)
     *
     * @return address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way).
     */
    virtual TissueCellPtr AddCell(TissueCellPtr pNewCell, const c_vector<double,DIM>& rCellDivisionVector, TissueCellPtr pParentCell=TissueCellPtr())=0;

    class Iterator; // Forward declaration; see below

    /**
     * Remove all cells labelled as dead.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @return number of cells removed
     */
    virtual unsigned RemoveDeadCells()=0;

    /**
     * Remove the Nodes (for cell-centre) or VertexElements (for cell-vertex) which
     * have been marked as deleted and update the correspondence with TissueCells.
     *
     * @param hasHadBirthsOrDeaths - a bool saying whether tissue has had Births Or Deaths
     */
    virtual void Update(bool hasHadBirthsOrDeaths=true)=0;

    /**
     * Find out how many cells of each mutation state there are
     *
     * @return The number of cells of each mutation state (evaluated at each visualizer output)
     * [0] = healthy count
     * [1] = APC one hit
     * [2] = APC two hit
     * [3] = beta catenin one hit
     */
    std::vector<unsigned> GetCellMutationStateCount();

    /**
     * Find out how many cells of each type there are
     *
     * @return The number of cells of each type (evaluated at each visualizer output)
     * [0] = STEM
     * [1] = TRANSIT
     * [2] = DIFFERENTIATED
     */
     const std::vector<unsigned>& rGetCellProliferativeTypeCount() const;

    /**
     * Find out how many cells in each cell cycle phase there are
     *
     * @return The number of cells of each phase (evaluated at each visualizer output)
     * [0] = G_ZERO_PHASE
     * [1] = G_ONE_PHASE
     * [2] = S_PHASE
     * [3] = G_TWO_PHASE
     * [4] = M_PHASE
     */
    const std::vector<unsigned>& rGetCellCyclePhaseCount() const;

    /**
     * Get the number of real cells.
     */
    unsigned GetNumRealCells();

    /**
     * Sets the Ancestor index of all the cells at this time to be the
     * same as their location index, can be used to trace clonal populations.
     */
    void SetCellAncestorsToLocationIndices();

    /**
     * Write cell ID data to mpCellIdFile.
     */
    void WriteCellIdDataToFile();

    /**
     * Loops over cells and makes a list of the ancestors that
     * are part of the tissue.
     *
     * @return remaining_ancestors  The size of this set tells you how many clonal populations remain.
     */
    std::set<unsigned> GetCellAncestors();

    /**
     * Get the cell corresponding to a given location index.
     *
     * Currently assumes there is one cell for each location index, and they are ordered identically in their vectors.
     * An assertion fails if not.
     *
     * @param index the location index
     *
     * @return the cell.
     */
    TissueCellPtr GetCellUsingLocationIndex(unsigned index);

    /**
     * Get the location index corresponding to a given cell.
     *
     * Currently assumes there is one cell for each location index, and they are ordered identically in their vectors.
     * An assertion fails if not.
     *
     * @param pCell the cell
     *
     * @return the location index.
     */
    unsigned GetLocationIndexUsingCell(TissueCellPtr pCell);

    /**
     * @return registry of cell properties used in this tissue.
     */
    boost::shared_ptr<CellPropertyRegistry> GetCellPropertyRegistry();

    /**
     * Set a default ordering on mutation states, so that existing tests don't need to
     * specify the old ordering explicitly.
     */
    void SetDefaultMutationStateOrdering();

    /**
     * Use an output file handler to create output files for visualizer and post-processing.
     *
     * @param rDirectory  pathname of the output directory, relative to where Chaste output is stored
     * @param cleanOutputDirectory  whether to delete the contents of the output directory prior to output file creation
     */
    virtual void CreateOutputFiles(const std::string& rDirectory,
                                   bool cleanOutputDirectory);

    /**
     * Write results from the current tissue state to output files.
     */
    virtual void WriteResultsToFiles();

    /**
     * Write the current time and node results to output files.
     */
    virtual void WriteTimeAndNodeResultsToFiles();

    /**
     * Calls GenerateCellResults() on each cell then calls WriteCellResultsToFiles().
     */
    virtual void GenerateCellResultsAndWriteToFiles()=0;

    /**
     * Generate results for a given cell in the current tissue state to output files.
     *
     * @param locationIndex location index of the cell
     * @param rCellProliferativeTypeCounter cell type counter
     * @param rCellCyclePhaseCounter cell cycle phase counter
     */
    virtual void GenerateCellResults(unsigned locationIndex,
                             std::vector<unsigned>& rCellProliferativeTypeCounter,
                             std::vector<unsigned>& rCellCyclePhaseCounter);

    /**
     * Write the current state of each cell to output files.

     * @param rCellProliferativeTypeCounter cell type counter
     * @param rCellCyclePhaseCounter cell cycle phase counter
     */
    void WriteCellResultsToFiles(std::vector<unsigned>& rCellProliferativeTypeCounter,
                                 std::vector<unsigned>& rCellCyclePhaseCounter);

    /**
     * Close any output files.
     */
    virtual void CloseOutputFiles();

    /**
     * @return mOutputCellIdData
     */
    bool GetOutputCellIdData();

    /**
     * @return mOutputCellMutationStates
     */
    bool GetOutputCellMutationStates();

    /**
     * @return mOutputCellAncestors
     */
    bool GetOutputCellAncestors();

    /**
     * @return mOutputCellProliferativeTypes
     */
    bool GetOutputCellProliferativeTypes();

    /**
     * @return mOutputCellVariables
     */
    bool GetOutputCellVariables();

    /**
     * @return mOutputCellCyclePhases
     */
    bool GetOutputCellCyclePhases();

    /**
     * @return mOutputCellAges
     */
    bool GetOutputCellAges();

    /**
     * @return mOutputCellVolumes
     */
    bool GetOutputCellVolumes();

    /**
     * Set mOutputCellIdData.
     * 
     * @param outputCellIdData the new value of mOutputCellIdData
     */
    void SetOutputCellIdData(bool outputCellIdData);

    /**
     * Set mOutputCellMutationStates.
     * 
     * @param outputCellMutationStates the new value of mOutputCellMutationStates
     */
    void SetOutputCellMutationStates(bool outputCellMutationStates);

    /**
     * Set mOutputCellAncestors.
     * 
     * @param outputCellAncestors the new value of mOutputCellAncestors
     */
    void SetOutputCellAncestors(bool outputCellAncestors);

    /**
     * Set mOutputCellProliferativeTypes.
     * 
     * @param outputCellProliferativeTypes the new value of mOutputCellProliferativeTypes
     */
    void SetOutputCellProliferativeTypes(bool outputCellProliferativeTypes);

    /**
     * Set mOutputCellVariables.
     * 
     * @param outputCellVariables the new value of mOutputCellVariables
     */
    void SetOutputCellVariables(bool outputCellVariables);

    /**
     * Set mOutputCellCyclePhases.
     * 
     * @param outputCellCyclePhases the new value of mOutputCellCyclePhases
     */
    void SetOutputCellCyclePhases(bool outputCellCyclePhases);

    /**
     * Set mOutputCellAges.
     * 
     * @param outputCellAges the new value of mOutputCellAges
     */
    void SetOutputCellAges(bool outputCellAges);

    /**
     * Set mOutputCellVolumes.
     * 
     * @param outputCellVolumes the new value of mOutputCellVolumes
     */
    void SetOutputCellVolumes(bool outputCellVolumes);

    /**
     * Iterator class allows one to iterate over cells in the tissue.
     * Dereferencing the iterator will give you the current cell.
     * There are also methods to get the node representing this cell,
     * and the location of that node.
     */
    class Iterator
    {
    public:

        /**
         * Dereference the iterator giving you a pointer to the current cell.
         */
        inline TissueCellPtr operator*();

        /**
         * Unusually for an iterator over a collection of pointers, this method
         * allows you to access the object pointed at, rather than the pointer itself.
         */
        inline TissueCellPtr operator->();

        /**
         * Comparison not-equal-to.
         *
         * @param rOther iterator with which comparison is made
         */
        inline bool operator!=(const AbstractTissue<DIM>::Iterator& rOther);

        /**
         * Prefix increment operator.
         */
        inline Iterator& operator++();

        /**
         * Constructor for a new iterator.
         *
         * @param rTissue the tissue
         * @param cellIter iterator over list of cells
         */
        Iterator(AbstractTissue& rTissue, std::list<TissueCellPtr>::iterator cellIter);

        /**
         * The iterator must have a virtual destructor.
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

        /** The tissue member. */
        AbstractTissue& mrTissue;

        /** Cell iterator member. */
        std::list<TissueCellPtr>::iterator mCellIter;
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

enum cell_colours
{
    STEM_COLOUR, // 0
    TRANSIT_COLOUR, // 1
    DIFFERENTIATED_COLOUR, // 2
    EARLY_CANCER_COLOUR, // 3
    LATE_CANCER_COLOUR, // 4
    LABELLED_COLOUR, // 5
    APOPTOSIS_COLOUR, // 6
    INVISIBLE_COLOUR // visualizer treats '7' as invisible
};

TEMPLATED_CLASS_IS_ABSTRACT_1_UNSIGNED(AbstractTissue)

//////////////////////////////////////////////////////////////////////////////
//         Iterator class implementation - most methods are inlined         //
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
TissueCellPtr AbstractTissue<DIM>::Iterator::operator*()
{
    assert(!IsAtEnd());
    return *mCellIter;
}

template<unsigned DIM>
TissueCellPtr AbstractTissue<DIM>::Iterator::operator->()
{
    assert(!IsAtEnd());
    return *mCellIter;
}

template<unsigned DIM>
bool AbstractTissue<DIM>::Iterator::operator!=(const AbstractTissue<DIM>::Iterator& rOther)
{
    return mCellIter != rOther.mCellIter;
}

template<unsigned DIM>
typename AbstractTissue<DIM>::Iterator& AbstractTissue<DIM>::Iterator::operator++()
{
    do
    {
        ++mCellIter;
    }
    while (!IsAtEnd() && !IsRealCell());

    return (*this);
}

template<unsigned DIM>
bool AbstractTissue<DIM>::Iterator::IsRealCell()
{
    return !( mrTissue.IsCellAssociatedWithADeletedLocation(*mCellIter) || (*this)->IsDead() );
}

template<unsigned DIM>
bool AbstractTissue<DIM>::Iterator::IsAtEnd()
{
    return mCellIter == mrTissue.rGetCells().end();
}

template<unsigned DIM>
AbstractTissue<DIM>::Iterator::Iterator(AbstractTissue& rTissue, std::list<TissueCellPtr>::iterator cellIter)
    : mrTissue(rTissue),
      mCellIter(cellIter)
{
    // The tissue can now return empty if it only has ghost nodes
    if (mrTissue.rGetCells().size() == 0)
    {
        mCellIter = mrTissue.rGetCells().end();
    }
    else
    {
        // Make sure we start at a real cell
        if (mCellIter == mrTissue.rGetCells().begin() && !IsRealCell())
        {
            ++(*this);
        }
    }
}

template<unsigned DIM>
typename AbstractTissue<DIM>::Iterator AbstractTissue<DIM>::Begin()
{
    return Iterator(*this, this->mCells.begin());
}

template<unsigned DIM>
typename AbstractTissue<DIM>::Iterator AbstractTissue<DIM>::End()
{
    return Iterator(*this, this->mCells.end());
}


#endif /*ABSTRACTTISSUE_HPP_*/
