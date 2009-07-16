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
#ifndef ABSTRACTTISSUE_HPP_
#define ABSTRACTTISSUE_HPP_

#include "TissueCell.hpp"
#include "OutputFileHandler.hpp"

#include <list>

#include <climits> // work around boost bug
#include <boost/serialization/access.hpp>
#include <boost/serialization/is_abstract.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/map.hpp>


/**
 * An abstract facade class encapsulating a tissue.
 *
 * Contains a group of cells and associated methods.
 *
 */
template<unsigned DIM>
class AbstractTissue
{
protected:

    /** List of cells */
    std::list<TissueCell> mCells;

    /** Map location (node or VertexElement) indices back to cells */
    std::map<unsigned, TissueCell*> mLocationCellMap;

    /** Map cells to location (node or VertexElement) indices */
    std::map<TissueCell*, unsigned> mCellLocationMap;
        
    /** Current cell mutation state counts */
    c_vector<unsigned, NUM_CELL_MUTATION_STATES> mCellMutationStateCount;

    /** Current cell type counts */
    c_vector<unsigned, NUM_CELL_TYPES> mCellTypeCount;

    /** Current cell cycle phase counts */
    c_vector<unsigned, 5> mCellCyclePhaseCount;

    /** Results file for node visualization */
    out_stream mpVizNodesFile;

    /** Results file for cell visualization */
    out_stream mpVizCellTypesFile;

    /** Results file for cell mutation states */
    out_stream mpCellMutationStatesFile;

    /** Results file for cell ancestors */
    out_stream mpCellAncestorsFile;

    /** Results file for cell types */
    out_stream mpCellTypesFile;

    /** Results file for cell cycle phases */
    out_stream mpCellCyclePhasesFile;

    /** Results file for cell variables */
    out_stream mpCellVariablesFile;

    /** Results file for cell ages */
    out_stream mpCellAgesFile;

    /** Results file for logged cell data. */
    out_stream mpCellIdFile;

    /** Whether the tissue contains a mesh */
    bool mTissueContainsMesh;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the facade.
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
        archive & mTissueContainsMesh;
    }

    /**
     * Check consistency of our internal data structures.
     */
    virtual void Validate()=0;

public:

    /**
     * Default constructor.
     *
     * @param rCells a vector of cells
     * @param locationIndices an optional vector of location indices that correspond to real cells
     */
    AbstractTissue(const std::vector<TissueCell>& rCells,
                   const std::vector<unsigned> locationIndices=std::vector<unsigned>());

    /**
     * Constructor for use by archiving - doesn't take in cells, since these are dealt
     * with by the serialize method.
     */
    AbstractTissue()
    {}

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
    std::list<TissueCell>& rGetCells();

    /**
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
     * @param pCell pointer to the cell
     *
     * @return the location of the cell
     */
    virtual c_vector<double, DIM> GetLocationOfCellCentre(TissueCell* pCell)=0;

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
     * @param rCell the cell
     *
     * @return whether a given cell is associated with a deleted node.
     */
    virtual bool IsCellAssociatedWithADeletedNode(TissueCell& rCell)=0;

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
     * @param rNewCell  the cell to add
     * @param newLocation  the position in space at which to put it
     * @param pParentCell pointer to a parent cell (if required)
     *
     * @return address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way).
     */
    virtual TissueCell* AddCell(TissueCell& rNewCell, c_vector<double,DIM> newLocation, TissueCell* pParentCell=NULL)=0;

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
     * [1] = labelled cells
     * [2] = APC one hit
     * [3] = APC two hit
     * [4] = beta catenin one hit
     */
    c_vector<unsigned, NUM_CELL_MUTATION_STATES> GetCellMutationStateCount();

     /**
     * Find out how many cells of each type there are
     *
     * @return The number of cells of each type (evaluated at each visualizer output)
     * [0] = STEM
     * [1] = TRANSIT
     * [2] = DIFFERENTIATED
     * [3] = APOPTOTIC
     */
    c_vector<unsigned, NUM_CELL_TYPES> GetCellTypeCount();

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
    c_vector<unsigned, 5> GetCellCyclePhaseCount();

    /**
     * Find if a given node is a ghost node. The abstract method always returns false
     * but is overridden in subclasses.
     *
     * @param index the global index of a specified node
     *
     * @return whether the node is a ghost node
     */
    virtual bool IsGhostNode(unsigned index);

    /**
     * Get the number of real cells.
     */
    unsigned GetNumRealCells();

    /**
     * Sets the Ancestor index of all the cells at this time to be the
     * same as their node index, can be used to trace clonal populations.
     */
    void SetCellAncestorsToNodeIndices();

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
    TissueCell& rGetCellUsingLocationIndex(unsigned index);

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
    unsigned GetLocationIndexUsingCell(TissueCell* pCell);

    /**
     * If the tissue contains a mesh, write this to file. For use by
     * the TissueSimulationArchiver. Must be overridden in each subclass
     * that contains a mesh.
     *
     * @param rArchiveDirectory directory in which archive is stored
     * @param rMeshFileName base name for mesh files
     */
    virtual void WriteMeshToFile(const std::string& rArchiveDirectory, const std::string& rMeshFileName);

    /**
     * Use an output file handler to create output files for visualizer and post-processing.
     *
     * @param rDirectory  pathname of the output directory, relative to where Chaste output is stored
     * @param rCleanOutputDirectory  whether to delete the contents of the output directory prior to output file creation
     * @param outputCellMutationStates  whether to create a cell mutation state results file
     * @param outputCellTypes  whether to create a cell type results file
     * @param outputCellVariables  whether to create a cell-cycle variable results file
     * @param outputCellCyclePhases  whether to create a cell-cycle phase results file
     * @param outputCellAncestors  whether to create a cell ancestor results file
     * @param outputCellAges whether to output cell age results
     */
    virtual void CreateOutputFiles(const std::string& rDirectory,
                                   bool rCleanOutputDirectory);

    /**
     * Write results from the current tissue state to output files.
     *
     * @param outputCellMutationStates  whether to output cell mutation state results
     * @param outputCellTypes  whether to output cell type results
     * @param outputCellVariables  whether to output cell-cycle variable results
     * @param outputCellCyclePhases  whether to output cell-cycle phase results
     * @param outputCellAncestors  whether to output cell ancestor results
     * @param outputCellAges whether to output cell age results
     */
    virtual void WriteResultsToFiles();

    /**
     * Write the current time and node results to output files.
     *
     * @param outputCellMutationStates  whether to output cell mutation state results
     * @param outputCellTypes  whether to output cell type results
     * @param outputCellVariables  whether to output cell-cycle variable results
     * @param outputCellCyclePhases  whether to output cell-cycle phase results
     * @param outputCellAncestors  whether to output cell ancestor results
     * @param outputCellAges whether to output cell age results
     * @param rCellTypeCounter cell type counter
     * @param rCellMutationStateCounter cell mutation state counter
     * @param rCellCyclePhaseCounter cell cycle phase counter
     */
    void WriteTimeAndNodeResultsToFiles(std::vector<unsigned>& rCellTypeCounter,
                                        std::vector<unsigned>& rCellMutationStateCounter,
                                        std::vector<unsigned>& rCellCyclePhaseCounter);

    /**
     * Generate results for a given cell in the current tissue state to output files.
     *
     * @param locationIndex location index of the cell
     * @param outputCellMutationStates  whether to output cell mutation state results
     * @param outputCellTypes  whether to output cell type results
     * @param outputCellVariables  whether to output cell-cycle variable results
     * @param outputCellCyclePhases  whether to output cell-cycle phase results
     * @param outputCellAncestors  whether to output cell ancestor results
     * @param outputCellAges whether to output cell age results
     * @param rCellTypeCounter cell type counter
     * @param rCellMutationStateCounter cell mutation state counter
     * @param rCellCyclePhaseCounter cell cycle phase counter
     */
    void GenerateCellResults(unsigned locationIndex,
    						 std::vector<unsigned>& rCellTypeCounter,
                             std::vector<unsigned>& rCellMutationStateCounter,
                             std::vector<unsigned>& rCellCyclePhaseCounter);

    /**
     * Write the current state of each cell to output files.
     *
     * @param outputCellMutationStates  whether to output cell mutation state results
     * @param outputCellTypes  whether to output cell type results
     * @param outputCellVariables  whether to output cell-cycle variable results
     * @param outputCellCyclePhases  whether to output cell-cycle phase results
     * @param outputCellAncestors  whether to output cell ancestor results
     * @param outputCellAges whether to output cell age results
     * @param rCellTypeCounter cell type counter
     * @param rCellMutationStateCounter cell mutation state counter
     * @param rCellCyclePhaseCounter cell cycle phase counter
     */
    void WriteCellResultsToFiles(std::vector<unsigned>& rCellTypeCounter,
                                 std::vector<unsigned>& rCellMutationStateCounter,
                                 std::vector<unsigned>& rCellCyclePhaseCounter);

    /**
     * Close any output files.
     *
     * @param outputCellMutationStates  whether a cell mutation state results file is open
     * @param outputCellTypes  whether a cell type results file is open
     * @param outputCellVariables  whether a cell-cycle variable results file is open
     * @param outputCellCyclePhases  whether a cell-cycle phase results file is open
     * @param outputCellAncestors  whether a cell ancestor results file is open
     * @param outputCellAges whether to output cell age results
     */
    virtual void CloseOutputFiles();

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
         * Dereference the iterator giving you a *reference* to the current cell.
         * Make sure to use a reference for the result to avoid copying cells unnecessarily.
         */
        inline TissueCell& operator*();

        /**
         * Member access from a pointer.
         */
        inline TissueCell* operator->();

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
        Iterator(AbstractTissue& rTissue, std::list<TissueCell>::iterator cellIter);

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
        std::list<TissueCell>::iterator mCellIter;
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
    INVISIBLE_COLOUR, // visualizer treats '7' as invisible
};

namespace boost
{
namespace serialization
{
/**
 * Since this abstract class is templated, we cannot use
 * the preprocessor macro BOOST_IS_ABSTRACT, and instead
 * must drop down to the underlying source code.
 */
template<unsigned DIM>
struct is_abstract<AbstractTissue<DIM> >
{
    /** The type that is an abstract class. */
    typedef mpl::bool_<true> type;
    /** The type is an abstract class, so value=true. */
    BOOST_STATIC_CONSTANT(bool, value=true);
};
}
}

//////////////////////////////////////////////////////////////////////////////
//         Iterator class implementation - most methods are inlined         //
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
TissueCell& AbstractTissue<DIM>::Iterator::operator*()
{
    assert(!IsAtEnd());
    return *mCellIter;
}

template<unsigned DIM>
TissueCell* AbstractTissue<DIM>::Iterator::operator->()
{
    assert(!IsAtEnd());
    return &(*mCellIter);
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
    return !( mrTissue.IsCellAssociatedWithADeletedNode(*mCellIter) || (*this)->IsDead() );
}

template<unsigned DIM>
bool AbstractTissue<DIM>::Iterator::IsAtEnd()
{
    return mCellIter == mrTissue.rGetCells().end();
}

template<unsigned DIM>
AbstractTissue<DIM>::Iterator::Iterator(AbstractTissue& rTissue, std::list<TissueCell>::iterator cellIter)
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
