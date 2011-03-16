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
#ifndef NODEBASEDCELLPOPULATION_HPP_
#define NODEBASEDCELLPOPULATION_HPP_

#include "AbstractCentreBasedCellPopulation.hpp"
#include "AbstractTetrahedralMesh.hpp" // for constructor which takes in a mesh
#include "BoxCollection.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>

/**
 * A NodeBasedCellPopulation is a CellPopulation consisting of only nodes in space with associated cells.
 * There are no elements and no mesh.
 */
template<unsigned DIM>
class NodeBasedCellPopulation : public AbstractCentreBasedCellPopulation<DIM>
{
    friend class TestNodeBasedCellPopulation;
    friend class TestBoxCollection;

protected:

    /** List of nodes. */
    std::vector<Node<DIM>* > mNodes;

    /** Indices of nodes that have been deleted, to be reused when adding new nodes. */
    std::vector<unsigned> mDeletedNodeIndices;

    /** Whether nodes have been added to the cell population. */
    bool mAddedNodes;

private:

    /** Pointer to a Node box collection */
    BoxCollection<DIM>* mpBoxCollection;

    /** Vector of minimal spatial positions in each dimension */
    c_vector<double, DIM> mMinSpatialPositions;

    /** Vector of maximal spatial positions in each dimension */
    c_vector<double, DIM> mMaxSpatialPositions;

    /** Node pairs for force calculations */
    std::set< std::pair<Node<DIM>*, Node<DIM>* > > mNodePairs;

    /**
     *  Whether to delete the nodes (taken in one of the constructors, defaults to true)
     */
    bool mDeleteNodes;

    /**
     * Mechanics cut off length.
     * Used in order to calculate the BoxCollection.
     */
    double mMechanicsCutOffLength;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * Note that serialization of the nodes is handled by load/save_construct_data,
     * so we don't actually have to do anything here except delegate to the base class.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCentreBasedCellPopulation<DIM> >(*this);

        Validate(); // paranoia
    }

    /**
     * Overridden AddNode() method.
     *
     * Add a new node to the cell population.
     *
     * @param pNewNode pointer to the new node
     * @return global index of new node in cell population
     */
    unsigned AddNode(Node<DIM>* pNewNode);

    /**
     * Move the node with a given index to a new point in space.
     *
     * @param nodeIndex the index of the node to be moved
     * @param rNewLocation the new target location of the node
     */
    void SetNode(unsigned nodeIndex, ChastePoint<DIM>& rNewLocation);

    /**
     * Check consistency of our internal data structures.
     */
    void Validate();

    /**
     * Method for Initially Splitting up cell population into neighbouring boxes, to decrease runtime.
     *
     * @param cutOffLength length of spring cut off between nodes
     * @param domainSize c_vector of size 2*dimension reads minX, maxX, minY, maxY, etc
     */
    void SplitUpIntoBoxes(double cutOffLength, c_vector<double, 2*DIM> domainSize);

    /**
     * Loops over nodes and sets mMinSpatialPositions and mMaxSpatialPositions
     */
    void FindMaxAndMin();

    /**
     * Overridden WriteVtkResultsToFile() method.
     */
    void WriteVtkResultsToFile();

public:

    /**
     * Default constructor.
     *
     * Note that the cell population will take responsibility for freeing the memory used by the nodes.
     *
     * @param nodes a vector of Nodes
     * @param rCells a vector of cells
     * @param locationIndices an optional vector of location indices that correspond to real cells
     * @param deleteNodes whether to delete nodes in destructor
     */
    NodeBasedCellPopulation(const std::vector<Node<DIM>* > nodes,
                    std::vector<CellPtr>& rCells,
                    const std::vector<unsigned> locationIndices=std::vector<unsigned>(),
                    bool deleteNodes=true);

    /**
     * Constructor for use by the archiving - doesn't take in cells, since these are
     * dealt with by the serialize method of our base class.
     *
     * Note that the cell population will take responsibility for freeing the memory used by the nodes.
     *
     * @param nodes a vector of Nodes
     * @param mechanicsCutOffLength the cut of length for mechanics uised to create the BoxedCollection to improve speed
     * @param deleteNodes whether to delete nodes in destructor
     */
    NodeBasedCellPopulation(const std::vector<Node<DIM>* > nodes, double mechanicsCutOffLength, bool deleteNodes=true);

    /**
     * Constructor which takes in a mesh and takes a copy of its nodes. The mesh is not
     * changed and no references to any of its data are created.
     *
     * This constructor is a helper constructor: it is generally easier for the user to
     * create a mesh than a set of nodes.
     *
     * @param rMesh any mesh.
     * @param rCells a vector of cells.
     */
    NodeBasedCellPopulation(const AbstractMesh<DIM,DIM>& rMesh,
                    std::vector<CellPtr>& rCells);

    /**
     * Destructor.
     *
     * Frees all our node memory.
     */
    ~NodeBasedCellPopulation();

    /**
     * @return the number of nodes in the cell population.
     */
    unsigned GetNumNodes();

    /**
     * Overridden GetNode() method.
     *
     * @param index global index of the specified node
     *
     * @return a pointer to the node with a given index.
     */
    Node<DIM>* GetNode(unsigned index);

    /**
     * Remove all cells labelled as dead.
     *
     * Note that after calling this method the cell population will be in an inconsistent state until
     * the equivalent of a 'remesh' is performed! So don't try iterating over cells or anything
     * like that.
     *
     * @return number of cells removed
     */
    unsigned RemoveDeadCells();

    /**
     * Reset the member variables #mDeletedNodeIndices, #mAddedNodes, #mNodePairs, and
     * #mpBoxCollection.
     */
    void Clear();

    /**
     * Remove nodes that have been marked as deleted and update the node cell map.
     *
     * @param hasHadBirthsOrDeaths whether cell population has had Births Or Deaths
     */
    void Update(bool hasHadBirthsOrDeaths=true);

    /**
     * Method for getting all nodes in the cell population.
     *
     * @return vector of Nodes
     */
    std::vector<Node<DIM>* >& rGetNodes();

    /**
     * Method for getting all nodes in the cell population (for archiving).
     *
     * @return vector of Nodes
     */
    const std::vector<Node<DIM>* >& rGetNodes() const;

    /**
     * @return pointer to a node box collection.
     */
    BoxCollection<DIM>* GetBoxCollection();

    /**
     * @return Node pairs for force calculation.
     */
    std::set< std::pair<Node<DIM>*, Node<DIM>* > >& rGetNodePairs();

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
     * @return mMechanicsCutOffLength
     */
    double GetMechanicsCutOffLength();

    /**
     * Set mMechanicsCutOffLength.
     * 
     * @param mechanicsCutOffLength  the new value of mMechanicsCutOffLength
     */
    void SetMechanicsCutOffLength(double mechanicsCutOffLength);

    /**
     * Overridden GetWidth() method.
     * 
     * Calculate the 'width' of any dimension of the cell population by computing
     * the maximum distance between any nodes in this dimension.
     *
     * @param rDimension a dimension (0,1 or 2)
     * @return The maximum distance between any nodes in this dimension.
     */
    double GetWidth(const unsigned& rDimension);

    /**
     * Overridden SetOutputCellVolumes() method.
     * 
     * Currently there is no facility for computing the volume associated
     * with each cell in a NodeBasedCellPopulation, so if this method is
     * called with outputCellVolumes = true, an exception is thrown.
     *
     * @param outputCellVolumes the new value of mOutputCellVolumes
     */
    void SetOutputCellVolumes(bool outputCellVolumes);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NodeBasedCellPopulation)

namespace boost {
namespace serialization {

/**
 * Non-intrusive serialization for Node - save method.
 */
template<class Archive, unsigned SPACE_DIM>
inline void save(
    Archive & ar,
    const Node<SPACE_DIM>& rNode,
    const unsigned int /* file_version */)
{
    // Save deleted flag
    const bool is_deleted = rNode.IsDeleted();
    ar << is_deleted;
}

/**
 * Non-intrusive serialization for Node - load method.
 */
template<class Archive, unsigned SPACE_DIM>
inline void load(
    Archive & ar,
    Node<SPACE_DIM>& rNode,
    const unsigned int /* file_version */)
{
    // Load deleted flag
    bool is_deleted;
    ar >> is_deleted;
    #define COVERAGE_IGNORE
    if (is_deleted)
    {
        rNode.MarkAsDeleted();
    }
    #undef COVERAGE_IGNORE
}

/**
 * Non-intrusive serialization for Node - serialize method.
 * This calls save or load as appropriate.
 */
template<class Archive, unsigned SPACE_DIM>
inline void serialize(
    Archive & ar,
    Node<SPACE_DIM>& rNode,
    const unsigned int file_version)
{
    boost::serialization::split_free(ar, rNode, file_version);
}

/**
 * Serialize information required to construct a Node.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const Node<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save the global index of the node
    const unsigned index = t->GetIndex();
    ar << index;

    // Save whether the node is a boundary node
    const bool is_boundary = t->IsBoundaryNode();
    ar << is_boundary;

    // Save the location of the node
    const c_vector<double, DIM>& r_loc = t->rGetLocation();
    for (unsigned i=0; i<DIM; i++)
    {
        ar << r_loc[i];
    }
}

/**
 * De-serialize constructor parameters and initialise a Node.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, Node<DIM> * t, const unsigned int file_version)
{
    // Load the global index of the node
    unsigned index;
    ar >> index;

    // Load whether the node is a boundary node
    bool is_boundary;
    ar >> is_boundary;

    // Load the location of the node
    c_vector<double, DIM> loc;
    for (unsigned i=0; i<DIM; i++)
    {
        ar >> loc[i];
    }

    // Invoke inplace constructor to initialise instance
    ::new(t)Node<DIM>(index, loc, is_boundary);
}

/**
 * Serialize information required to construct a NodeBasedCellPopulation.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const NodeBasedCellPopulation<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    ar & t->rGetNodes();
//    const double cut_off_length = t->GetMechanicsCutOffLength();
//    ar << cut_off_length;

}

/**
 * De-serialize constructor parameters and initialise a NodeBasedCellPopulation.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, NodeBasedCellPopulation<DIM> * t, const unsigned int file_version)
{
    // Load the nodes
    std::vector<Node<DIM>* > nodes;
    ar >> nodes;

    // load the cut of distance
    double cut_off_length = 1.0; //\todo this is temporary till we fix archiving
//    ar >> cut_off_length;

    // Invoke inplace constructor to initialize instance
    ::new(t)NodeBasedCellPopulation<DIM>(nodes, cut_off_length);
}

}} // close namespaces

#endif /*NODEBASEDCELLPOPULATION_HPP_*/
