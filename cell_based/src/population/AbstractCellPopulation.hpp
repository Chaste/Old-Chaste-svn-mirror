/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef ABSTRACTCELLPOPULATION_HPP_
#define ABSTRACTCELLPOPULATION_HPP_

#include "Cell.hpp"
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
#include "Identifiable.hpp"

// Needed here to avoid serialization errors (on Boost<1.37)
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "ApoptoticCellProperty.hpp"
#include "CellLabel.hpp"

/**
 * An abstract facade class encapsulating a cell population.
 *
 * Contains a group of cells and associated methods.
 *
 */
template<unsigned DIM>
class AbstractCellPopulation : public Identifiable
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
    std::list<CellPtr> mCells;

    /** Map location (node or VertexElement) indices back to cells */
    std::map<unsigned, CellPtr> mLocationCellMap;

    /** Map cells to location (node or VertexElement) indices */
    std::map<Cell*, unsigned> mCellLocationMap;

    /** Current cell type counts */
    std::vector<unsigned> mCellProliferativeTypeCount;

    /** Current cell cycle phase counts */
    std::vector<unsigned> mCellCyclePhaseCount;

    /** Population centroid */
    c_vector<double, DIM> mCentroid;

    /** Results file for node visualization */
    out_stream mpVizNodesFile;

    /** Results file for cell visualization */
    out_stream mpVizCellProliferativeTypesFile;

    /** Results file for cell visualization */
    out_stream mpVizCellProliferativePhasesFile;

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
        
    /** Results file for cell volume (in 3D) or area (in 2D) data. */
    out_stream mpCellVolumesFile;

    /** Results file for boundary nodes. */
    out_stream mpVizBoundaryNodesFile;

    /** A cache of where the results are going (used for VTK writer). */
    std::string mDirPath;

    /** Meta results file for VTK. */
    out_stream mpVtkMetaFile;

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
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     */
    virtual void Validate()=0;

    /**
     * Write the current results to mpVtkMetaFile.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     */
    virtual void WriteVtkResultsToFile()=0;

    /**
     * Constructor for use by archiving only. Please use the other constructor.
     *
     * Doesn't take in cells, since these are dealt with by the serialize method.
     */
    AbstractCellPopulation()
    {
    }

public:

    /**
     * AbstractCellPopulation Constructor.
     *
     * @note Warning: the passed-in vector of cells will be emptied, even if the constructor
     * throws an exception!
     *
     * @param rCells a vector of cells.  Copies of the cells will be stored in the cell population,
     *     and the passed-in vector cleared.
     * @param locationIndices an optional vector of location indices that correspond to real cells
     */
    AbstractCellPopulation(std::vector<CellPtr>& rCells,
                           const std::vector<unsigned> locationIndices=std::vector<unsigned>());

    /**
     * Base class with virtual methods needs a virtual destructor.
     */
    virtual ~AbstractCellPopulation();

    /**
     * Initialise each cell's cell-cycle model.
     */
    void InitialiseCells();

    /**
     * @return reference to mCells.
     */
    std::list<CellPtr>& rGetCells();

    /**
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @return the number of nodes in the cell population.
     */
    virtual unsigned GetNumNodes()=0;

    /**
     * Find where a given cell is in space.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param pCell the cell
     * @return the location of the cell
     */
    virtual c_vector<double, DIM> GetLocationOfCellCentre(CellPtr pCell)=0;

    /**
     * Get a pointer to the node with a given index.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param index  global index of the specified node
     * @return a pointer to the node with a given index.
     */
    virtual Node<DIM>* GetNode(unsigned index)=0;

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
     * @return whether a given cell is associated with a deleted
     *         node (cell-centre models) or element (vertex models).
     */
    virtual bool IsCellAssociatedWithADeletedLocation(CellPtr pCell)=0;

    /**
     * Add a new cell to the cell population.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param pNewCell  the cell to add
     * @param rCellDivisionVector  a vector providing information regarding how the cell division should occur
     *     (for cell-centre cell populations, this vector is the position of the daughter cell; for vertex cell populations it
     *      can be used by any subclass of CellBasedSimulation to as a means of dictating the axis along which
     *      the parent cell divides)
     * @param pParentCell pointer to a parent cell (if required)
     *
     * @return address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way).
     */
    virtual CellPtr AddCell(CellPtr pNewCell, const c_vector<double,DIM>& rCellDivisionVector, CellPtr pParentCell=CellPtr())=0;

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
     * have been marked as deleted and update the correspondence with Cells.
     *
     * @param hasHadBirthsOrDeaths - a bool saying whether cell population has had Births Or Deaths
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
     * are part of the cell population.
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
    CellPtr GetCellUsingLocationIndex(unsigned index);

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
    unsigned GetLocationIndexUsingCell(CellPtr pCell);

    /**
     * @return registry of cell properties used in this cell population.
     */
    boost::shared_ptr<CellPropertyRegistry> GetCellPropertyRegistry();

    /**
     * Set a default ordering on mutation states, so that existing tests don't need to
     * specify the old ordering explicitly.
     */
    void SetDefaultMutationStateOrdering();

    /**
     * Calculate the 'width' of any dimension of the cell population.
     *      
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rDimension a dimension (0,1 or 2)
     * @return The maximum distance between any nodes in this dimension.
     */
    virtual double GetWidth(const unsigned& rDimension)=0;

    /**
     * Given a node index, returns the set of neighbouring node indices.
     *      
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param index the node index
     * @return the set of neighbouring node indices.
     */
    virtual std::set<unsigned> GetNeighbouringNodeIndices(unsigned index)=0;

    /**
     * Returns the centroid of the cell population.
     */
    c_vector<double, DIM> GetCentroidOfCellPopulation();

    /**
     * Use an output file handler to create output files for visualizer and post-processing.
     *
     * @param rDirectory  pathname of the output directory, relative to where Chaste output is stored
     * @param cleanOutputDirectory  whether to delete the contents of the output directory prior to output file creation
     */
    virtual void CreateOutputFiles(const std::string& rDirectory,
                                   bool cleanOutputDirectory);

    /**
     * Write results from the current cell population state to output files.
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
     * Generate results for a given cell in the current cell population state to output files.
     *
     * @param locationIndex location index of the cell
     * @param rCellProliferativeTypeCounter cell type counter
     * @param rCellCyclePhaseCounter cell cycle phase counter
     */
    virtual void GenerateCellResults(unsigned locationIndex,
                             std::vector<unsigned>& rCellProliferativeTypeCounter,
                             std::vector<unsigned>& rCellCyclePhaseCounter);

    /**
     * Write the current volume of each cell to file.      
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     */
    virtual void WriteCellVolumeResultsToFile()=0;

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
     * Outputs CellPopulation used in the simulation to file and then calls OutputCellPopulationParameters to output all relevant parameters.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationInfo(out_stream& rParamsFile);

    /**
     * Outputs CellPopulation parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellPopulationParameters(out_stream& rParamsFile)=0;

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
     * Iterator class allows one to iterate over cells in the cell population.
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
        inline CellPtr operator*();

        /**
         * Unusually for an iterator over a collection of pointers, this method
         * allows you to access the object pointed at, rather than the pointer itself.
         */
        inline CellPtr operator->();

        /**
         * Comparison not-equal-to.
         *
         * @param rOther iterator with which comparison is made
         */
        inline bool operator!=(const AbstractCellPopulation<DIM>::Iterator& rOther);

        /**
         * Prefix increment operator.
         */
        inline Iterator& operator++();

        /**
         * Constructor for a new iterator.
         *
         * @param rCellPopulation the cell population
         * @param cellIter iterator over list of cells
         */
        Iterator(AbstractCellPopulation& rCellPopulation, std::list<CellPtr>::iterator cellIter);

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

        /** The cell population member. */
        AbstractCellPopulation& mrCellPopulation;

        /** Cell iterator member. */
        std::list<CellPtr>::iterator mCellIter;
    };

    /**
     * @return iterator pointing to the first cell in the cell population
     */
    Iterator Begin();

    /**
     * @return iterator pointing to one past the last cell in the cell population
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

TEMPLATED_CLASS_IS_ABSTRACT_1_UNSIGNED(AbstractCellPopulation)

//////////////////////////////////////////////////////////////////////////////
//         Iterator class implementation - most methods are inlined         //
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
CellPtr AbstractCellPopulation<DIM>::Iterator::operator*()
{
    assert(!IsAtEnd());
    return *mCellIter;
}

template<unsigned DIM>
CellPtr AbstractCellPopulation<DIM>::Iterator::operator->()
{
    assert(!IsAtEnd());
    return *mCellIter;
}

template<unsigned DIM>
bool AbstractCellPopulation<DIM>::Iterator::operator!=(const AbstractCellPopulation<DIM>::Iterator& rOther)
{
    return mCellIter != rOther.mCellIter;
}

template<unsigned DIM>
typename AbstractCellPopulation<DIM>::Iterator& AbstractCellPopulation<DIM>::Iterator::operator++()
{
    do
    {
        ++mCellIter;
    }
    while (!IsAtEnd() && !IsRealCell());

    return (*this);
}

template<unsigned DIM>
bool AbstractCellPopulation<DIM>::Iterator::IsRealCell()
{
    return !( mrCellPopulation.IsCellAssociatedWithADeletedLocation(*mCellIter) || (*this)->IsDead() );
}

template<unsigned DIM>
bool AbstractCellPopulation<DIM>::Iterator::IsAtEnd()
{
    return mCellIter == mrCellPopulation.rGetCells().end();
}

template<unsigned DIM>
AbstractCellPopulation<DIM>::Iterator::Iterator(AbstractCellPopulation& rCellPopulation, std::list<CellPtr>::iterator cellIter)
    : mrCellPopulation(rCellPopulation),
      mCellIter(cellIter)
{
    // The cell population can now return empty if it only has ghost nodes
    if (mrCellPopulation.rGetCells().empty())
    {
        mCellIter = mrCellPopulation.rGetCells().end();
    }
    else
    {
        // Make sure we start at a real cell
        if (mCellIter == mrCellPopulation.rGetCells().begin() && !IsRealCell())
        {
            ++(*this);
        }
    }
}

template<unsigned DIM>
typename AbstractCellPopulation<DIM>::Iterator AbstractCellPopulation<DIM>::Begin()
{
    return Iterator(*this, this->mCells.begin());
}

template<unsigned DIM>
typename AbstractCellPopulation<DIM>::Iterator AbstractCellPopulation<DIM>::End()
{
    return Iterator(*this, this->mCells.end());
}


#endif /*ABSTRACTCELLPOPULATION_HPP_*/
