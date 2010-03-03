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
#ifndef NODEBASEDTISSUE_HPP_
#define NODEBASEDTISSUE_HPP_

#include "AbstractCellCentreBasedTissue.hpp"
#include "AbstractTetrahedralMesh.hpp" // for constructor which takes in a mesh
#include "BoxCollection.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>

/**
 * A NodeBasedTissue is a Tissue consisting of only nodes in space with associated cells.
 * There are no elements and no mesh.
 */
template<unsigned DIM>
class NodeBasedTissue : public AbstractCellCentreBasedTissue<DIM>
{
    friend class TestNodeBasedTissue;
    friend class TestBoxCollection;

protected:

    /** List of nodes. */
    std::vector<Node<DIM>* > mNodes;

    /** Indices of nodes that have been deleted, to be reused when adding new nodes. */
    std::vector<unsigned> mDeletedNodeIndices;

    /** Whether nodes have been added to the tissue. */
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
        archive & boost::serialization::base_object<AbstractCellCentreBasedTissue<DIM> >(*this);

        Validate(); // paranoia
    }

    /**
     * Overridden AddNode() method.
     *
     * Add a new node to the tissue.
     *
     * @param pNewNode pointer to the new node
     * @return global index of new node in tissue
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
     * Method for Initially Splitting up tissue into neighbouring boxes, to decrease runtime.
     *
     * @param cutOffLength length of spring cut off between nodes
     * @param domainSize c_vector of size 2*dimension reads minX, maxX, minY, maxY, etc
     */
    void SplitUpIntoBoxes(double cutOffLength, c_vector<double, 2*DIM> domainSize);

    /**
     * Loops over nodes and sets mMinSpatialPositions and mMaxSpatialPositions
     */
    void FindMaxAndMin();

public:

    /**
     * Default constructor.
     *
     * Note that the tissue will take responsibility for freeing the memory used by the nodes.
     *
     * @param nodes a vector of Nodes
     * @param rCells a vector of TissueCells
     * @param locationIndices an optional vector of location indices that correspond to real cells
     * @param deleteNodes whether to delete nodes in destructor
     */
    NodeBasedTissue(const std::vector<Node<DIM>* > nodes,
                    const std::vector<TissueCell>& rCells,
                    const std::vector<unsigned> locationIndices=std::vector<unsigned>(),
                    bool deleteNodes=true);

    /**
     * Constructor for use by the archiving - doesn't take in cells, since these are
     * dealt with by the serialize method of our base class.
     *
     * Note that the tissue will take responsibility for freeing the memory used by the nodes.
     *
     * @param nodes a vector of Nodes
     * @param deleteNodes whether to delete nodes in destructor
     */
    NodeBasedTissue(const std::vector<Node<DIM>* > nodes, bool deleteNodes=true);

    /**
     * Constructor which takes in a mesh and takes a copy of its nodes. The mesh is not
     * changed and no references to any of its data are created.
     *
     * This constructor is a helper constructor: it is generally easier for the user to
     * create a mesh than a set of nodes.
     *
     * @param rMesh any mesh.
     * @param rCells a vector of TissueCells.
     */
    NodeBasedTissue(const AbstractMesh<DIM,DIM>& rMesh,
                    const std::vector<TissueCell>& rCells);

    /**
     * Destructor.
     *
     * Frees all our node memory.
     */
    ~NodeBasedTissue();

    /**
     * @return the number of nodes in the tissue.
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
     * Note that after calling this method the tissue will be in an inconsistent state until
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
     * @param hasHadBirthsOrDeaths whether tissue has had Births Or Deaths
     */
    void Update(bool hasHadBirthsOrDeaths=true);

    /**
     * Method for getting all nodes in the tissue.
     *
     * @return vector of Nodes
     */
    std::vector<Node<DIM>* >& rGetNodes();

    /**
     * Method for getting all nodes in the tissue (for archiving).
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
};


#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NodeBasedTissue)

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
 * Serialize information required to construct a NodeBasedTissue.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const NodeBasedTissue<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    ar & t->rGetNodes();
}

/**
 * De-serialize constructor parameters and initialise a NodeBasedTissue.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, NodeBasedTissue<DIM> * t, const unsigned int file_version)
{
    // Load the nodes
    std::vector<Node<DIM>* > nodes;
    ar >> nodes;

    // Invoke inplace constructor to initialize instance
    ::new(t)NodeBasedTissue<DIM>(nodes);
}

}} // close namespaces


#endif /*NODEBASEDTISSUE_HPP_*/
