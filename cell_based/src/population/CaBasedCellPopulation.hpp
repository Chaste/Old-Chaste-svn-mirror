/*

Copyright (c) 2005-2012, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef CABASEDCELLPOPULATION_HPP_
#define CABASEDCELLPOPULATION_HPP_

#include <vector>
#include <set>
#include <string>

#include "UblasVectorInclude.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

// Needed here to avoid serialization errors (on Boost<1.37)
#include "WildTypeCellMutationState.hpp"

#include "AbstractOnLatticeCellPopulation.hpp"
#include "TetrahedralMesh.hpp"
#include "AbstractCaUpdateRule.hpp"

template<unsigned DIM>
class AbstractCaUpdateRule; // Circular definition

/**
 * A facade class encapsulating a lattice-based cell population.
 *
 * Contains a group of cells and maintains the associations between cells and
 * nodes in the mesh which correspond to lattice sites.
 */
template<unsigned DIM>
class CaBasedCellPopulation : public AbstractOnLatticeCellPopulation<DIM>
{
    friend class TestCaBasedCellPopulation;

private:

    /** Reference to the mesh associated with the cell population. */
    TetrahedralMesh<DIM, DIM>& mrMesh;

    /** The update rules used to determine the new location of the cells. */
    std::vector<boost::shared_ptr<AbstractCaUpdateRule<DIM> > > mUpdateRuleCollection;

    /** Records whether a node is an empty site or not. */
    std::vector<bool> mIsEmptySite;

    /**
     * Whether to only search the next nearest neighbours for an empty site when a cell
     * is going to divide.
     * Defaults to false in the constructor.
     */
    bool mOnlyUseNearestNeighboursForDivision;

    /**
     * Whether to implement von Neumann neighbourhoods for dividing and moving cells.
     * These neighbourhoods only correspond to the N, S, E, W neighbours.
     * Defaults to false in the constructor.
     */
    bool mUseVonNeumannNeighbourhoods;

    /**
     * Set the empty sites by taking in a set of which nodes indices are empty sites.
     *
     * @param rEmptySiteIndices set of node indices corresponding to empty sites
     */
    void SetEmptySites(const std::set<unsigned>& rEmptySiteIndices);

    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * Note that serialization of the mesh and cells is handled by load/save_construct_data.
     *
     * Note also that member data related to writers is not saved - output must
     * be set up again by the caller after a restart.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOnLatticeCellPopulation<DIM> >(*this);
        archive & mUpdateRuleCollection;
        archive & mIsEmptySite;
        archive & mOnlyUseNearestNeighboursForDivision;
        archive & mUseVonNeumannNeighbourhoods;
    }

    /**
     * Check consistency of our internal data structures. Each node must
     * have a cell associated with it.
     */
    void Validate();

    /**
     * Overridden WriteVtkResultsToFile() method.
     */
    void WriteVtkResultsToFile();

public:

    /**
     * Create a new cell population facade from a mesh and collection of cells.
     *
     * There must be precisely 1 cell for each node of the mesh.
     *
     * @param rMesh a mutable tetrahedral mesh
     * @param rCells CellPtrs corresponding to the nodes of the mesh
     * @param locationIndices an optional vector of location indices that correspond to real cells
     * @param deleteMesh set to true if you want the cell population to free the mesh memory on destruction
     *            (defaults to false)
     * @param validate whether to validate the cell population (defaults to false)
     */
    CaBasedCellPopulation(TetrahedralMesh<DIM, DIM>& rMesh,
                       std::vector<CellPtr>& rCells,
                       const std::vector<unsigned> locationIndices=std::vector<unsigned>(),
                       bool deleteMesh=false,
                       bool validate=true);

    /**
     * Constructor for use by the de-serializer.
     *
     * @param rMesh a tetrahedral mesh.
     */
    CaBasedCellPopulation(TetrahedralMesh<DIM, DIM>& rMesh);

    /**
     * Destructor.
     */
    ~CaBasedCellPopulation();

    /**
     * @return reference to  mrMesh.
     */
    TetrahedralMesh<DIM, DIM>& rGetMesh();

    /**
     * @return const reference to mrMesh (used in archiving).
     */
    const TetrahedralMesh<DIM, DIM>& rGetMesh() const;

    /**
     * Overridden UpdateCellLocations() method.
     *
     * @param dt time step
     */
    void UpdateCellLocations(double dt);

    /**
     * Set mOnlyUseNearestNeighboursForDivision.
     *
     * @param onlyUseNearestNeighboursForDivision whether to only search the next nearest neighbours for an
     *     empty site when a cell is going to divide.
     */
    void SetOnlyUseNearestNeighboursForDivision(bool onlyUseNearestNeighboursForDivision);

    /**
     * Get mOnlyUseNearestNeighboursForDivision.
     */
    bool GetOnlyUseNearestNeighboursForDivision();

    /**
     * Add an update rule to be used in this simulation.
     *
     * @param pUpdateRule pointer to an update rule
     */
    void AddUpdateRule(boost::shared_ptr<AbstractCaUpdateRule<DIM> > pUpdateRule);

    /**
     * Get the collection of update rules to be used in this simulation.
     *
     * @return the update rule collection
     */
    const std::vector<boost::shared_ptr<AbstractCaUpdateRule<DIM> > >& rGetUpdateRuleCollection() const;

    /**
     * Set method for mUseVonNeumannNeighbourhoods.
     *
     * @param useVonNeumannNeighbourhoods whether to use von Neumann neighbourhoods
     */
    void SetUseVonNeumannNeighbourhoods(bool useVonNeumannNeighbourhoods);

    /**
     * @return mUseVonNeumannNeighbourhoods
     */
    bool GetUseVonNeumannNeighbourhoods();

    /**
     * @return mIsEmptySite.
     */
    std::vector<bool>& rGetEmptySites();

    /**
     * Find if a given node is an empty site. The abstract method always returns false
     * but is overridden in subclasses.
     *
     * @param index the global index of a specified node
     *
     * @return whether the node is an empty site
     */
    bool IsEmptySite(unsigned index);

    /**
     * @return the indices of those nodes that are empty sites.
     */
    std::set<unsigned> GetEmptySiteIndices();

    /**
     * Overridden GetNumNodes() method.
     *
     * @return the number of nodes in the cell population.
     */
    unsigned GetNumNodes();

    /**
     * Overridden GetNode() method.
     *
     * @param index  global index of the specified Node
     *
     * @return pointer to the Node with given index.
     */
    virtual Node<DIM>* GetNode(unsigned index);

    /**
     * Overridden AddCell() method.
     *
     * Add a new cell to the cell population and update mIsEmptySite.
     *
     * This method finds the nearest empty site to move to in the admissable directions (see diagram below).
     *
     * The bool mOnlyUseNearestNeighboursForDivision is defined in the cell population constructor and
     * says whether to only search the immediate nearest neighbours when dividing.
     *
     * If there is more than one empty site at the same distance (i.e. degree) from the parent cell,
     * then the method chooses randomly between them.
     *
     * If for a given degree there are not any empty sites, then it checks for empty sites at degree+1 (provided
     * mOnlyUseNearestNeighboursForDivision = false.
     *
     * This diagram indicates the possible directions and degrees of neighbouring sites from parent cell X:
     *
     *            3 o o 3 o o 3
     *            o 2 o 2 o 2 o
     *            o o 1 1 1 o o
     *            3 2 1 X 1 2 3
     *            o o 1 1 1 o o
     *            o 2 o 2 o 2 o
     *            3 o o 3 o o 3
     *
     * Once it has located an empty site, neighbouring cells are pushed along in that direction until the 1st degree neighbouring
     * site is free, and the new daughter cell is put there.
     *
     * If there are no free sites in any direction, currently the code throws an exception.
     *
     * @param pNewCell  the cell to add
     * @param rCellDivisionVector  the position in space at which to put it
     * @param pParentCell pointer to a parent cell, if required (defaults to NULL)
     *
     * @return address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way)
     */
    virtual CellPtr AddCell(CellPtr pNewCell, const c_vector<double,DIM>& rCellDivisionVector, CellPtr pParentCell=CellPtr());

    /**
     * Same method as GetNeighbouringNodeIndices() above, but returns an ORDERED vector
     * of neighbouring nodes, in the order N, NW, W, SW, S, SE, E, NE. This is needed
     * to make GetNthDegreeNeighbouringNodeIndices work, as this exploits the ordering
     * of the neighbours.
     *
     * @param nodeIndex global index of the node of interest
     */
    std::vector<unsigned> GetNeighbouringNodeIndicesVector(unsigned nodeIndex);

    /**
     * Locate the sites neighbouring a site (this version is a Moore neighbourhood)
     * which are free (i.e. are empty sites).
     *
     * Note: This dictates the geometry of the cell population and the type of neighbourhood
     * used and can be overridden to use different neighbourhoods or geometries.
     *
     * @param nodeIndex global index of the node of interest
     *
     * @return set of available nodes
     */
    std::set<unsigned> GetFreeNeighbouringNodeIndices(unsigned nodeIndex);

    /**
     * Return the maximum degree that is permissable for a given node in each direction
     * for which this maximum degree is non-zero.
     *
     * The vector uses the ordering N, NW, W, SW, S, SE, E, NE (i.e. anticlockwise).
     *
     * This degree corresponds to the number of nodes from the given node up to and
     * including the boundary node in each direction.
     *
     * @param nodeIndex index of the node of interest
     */
    std::vector<unsigned> GetMaximumDegreeInEachDirection(unsigned nodeIndex);

    /**
     * Locate the sites in n-th degree neighbouring sites (this version is a Moore
     * neighbourhood). Note: This dictates the geometry of the cell population and the type
     * of neighbourhood used and can be overridden to use cell population neighbourhoods
     * or geometries.
     *
     * @param nodeIndex global index of the node of interest
     * @param degree how far away to look for neighbours, i.e. 1 would be the nearest neighbours
     *
     * @return set of available nodes
     */
    std::set<unsigned> GetNthDegreeNeighbouringNodeIndices(unsigned nodeIndex, unsigned degree);

    /**
     * Remove all cells that are labelled as dead.
     *
     * Note that this now calls MutableMesh::DeleteNodePriorToReMesh()
     * and therefore a ReMesh(map) must be called before any element
     * information is used.
     *
     * Note also that after calling this method the cell population will be in an inconsistent state until
     * Update() is called! So don't try iterating over cells or anything like that.
     *
     * @return number of cells removed.
     */
    unsigned RemoveDeadCells();

    /**
     * Overridden Update(bool hasHadBirthsOrDeaths) method.
     * Fixes up the mappings between cells and nodes.
     *
     * @param hasHadBirthsOrDeaths whether there has been any cell division or cell
     *                             death prior to the update (defaults to true)
     */
    void Update(bool hasHadBirthsOrDeaths=true);

    /**
     * Overridden WriteCellVolumeResultsToFile() method.
     */
    void WriteCellVolumeResultsToFile();

    /**
     * Overridden GetVolumeOfCell() method.
     * 
     * @param pCell boost shared pointer to a cell
     */
    double GetVolumeOfCell(CellPtr pCell);

    /**
     * Overridden GenerateCellResults() method.
     * Generate results for a given cell in the current cell population state to output files.
     *
     * @param locationIndex location index of the cell
     * @param rCellProliferativeTypeCounter cell type counter
     * @param rCellCyclePhaseCounter cell cycle phase counter
     */
    void GenerateCellResults(unsigned locationIndex,
                             std::vector<unsigned>& rCellProliferativeTypeCounter,
                             std::vector<unsigned>& rCellCyclePhaseCounter);

    /**
     * Calls GenerateCellResults() on each cell then calls WriteCellResultsToFiles().
     */
    virtual void GenerateCellResultsAndWriteToFiles();

    /**
     * Helper method for establishing if a cell is real.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * Here we don't allow deleted nodes so this always returns false.
     *
     * @param pCell the cell
     *
     * @return whether a given cell is associated with a deleted
     *         node (cell-centre models) or element (vertex models).
     */
    bool IsCellAssociatedWithADeletedLocation(CellPtr pCell);

    /**
     * Method to update the location of a particular cell
     *
     * @param pCell the cell to move
     * @param newLocationIndex the location to move the cell to
     */
    void MoveCell(CellPtr pCell, unsigned newLocationIndex);

    /**
     * Find where a given cell is in space.
     *
     * @param pCell the cell
     *
     * @return the location of the cell
     */
    c_vector<double, DIM> GetLocationOfCellCentre(CellPtr pCell);

    /**
     * Outputs CellPopulation parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationParameters(out_stream& rParamsFile);

    /**
     * Overridden GetWidth() method.
     *
     * Calculate the 'width' of any dimension of the cell population by calling
     * GetWidth() on the mesh.
     *
     * @param rDimension a dimension (0,1 or 2)
     * @return The maximum distance between any nodes in this dimension.
     */
    double GetWidth(const unsigned& rDimension);

    /**
     * Overridden GetNeighbouringNodeIndices() method.
     *
     * Locate the sites neighbouring a site (this version is a Moore neighbourhood).
     * Note: This dictates the geometry of the cell population and the type of neighbourhood
     * used and can be overridden to use different neighbourhoods or geometries.
     *
     * @param index the node index
     * @return the set of neighbouring node indices.
     */
    std::set<unsigned> GetNeighbouringNodeIndices(unsigned index);
};
#undef COVERAGE_IGNORE // Avoid prototypes being treated as code by gcov

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CaBasedCellPopulation)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a CaBasedCellPopulation.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const CaBasedCellPopulation<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const TetrahedralMesh<DIM,DIM>* p_mesh = &(t->rGetMesh());
    ar & p_mesh;
}

/**
 * De-serialize constructor parameters and initialise a CellPopulation.
 * Loads the mesh from separate files.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, CaBasedCellPopulation<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    TetrahedralMesh<DIM,DIM>* p_mesh;
    ar >> p_mesh;

    // Invoke inplace constructor to initialise instance
    ::new(t)CaBasedCellPopulation<DIM>(*p_mesh);
}
}
} // namespace ...

#endif /*CABASEDCELLPOPULATION_HPP_*/
