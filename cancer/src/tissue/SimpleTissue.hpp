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
#ifndef SIMPLETISSUE_HPP_
#define SIMPLETISSUE_HPP_

#include "AbstractTissue.hpp"
#include "ConformingTetrahedralMesh.hpp"

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>


template<unsigned DIM>
class SimpleTissue : public AbstractTissue<DIM>
{
    friend class TestSimpleTissue;
private:

    /** List of nodes */
    std::vector<Node<DIM> > mNodes;    
    
    friend class boost::serialization::access;
    /**
     * Serialize the facade.
     * 
     * Note that serialization of the nodes is handled by load/save_construct_data,
     * so we don't actually have to do anything here except delegate to the base class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractTissue<DIM> >(*this);
        
        Validate(); // paranoia
    }
            
    /** 
     * Add a new node to the tissue. 
     */
    unsigned AddNode(Node<DIM> *pNewNode);
    
    /** 
     * Move the node with a given index to a new point in space.
     */
    void SetNode(unsigned index, ChastePoint<DIM> point);
    

public:

    SimpleTissue(const std::vector<Node<DIM> >& rNodes, const std::vector<TissueCell>& rCells);

    /**
     * Constructor for use by the archiving - doesn't take in cells, since these are
     * dealt with by the serialize method of our base class.
     */
    SimpleTissue(const std::vector<Node<DIM> >& rNodes);
    
    ~SimpleTissue() 
    {}

    /** 
     * Get the number of nodes in the tissue.
     */
    unsigned GetNumNodes();
    
    /**
     * Get a pointer to the node with a given index.
     */  
    Node<DIM>* GetNode(unsigned index);
    
    /**
     * Add a new cell to the tissue.
     * @param cell  the cell to add
     * @param newLocation  the position in space at which to put it
     * @returns address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way)
     */
    TissueCell* AddCell(TissueCell cell, c_vector<double,DIM> newLocation);
    
    /** 
     * Remove all cells labelled as dead. 
     * 
     * Note that after calling this method the tissue will be in an inconsistent state until 
     * the equivalent of a 'remesh' is performed! So don't try iterating over cells or anything 
     * like that.
     * 
     *  @return number of cells removed
     */
    unsigned RemoveDeadCells();
    
    /** 
     * Remove nodes that have been marked as deleted and update the node cell map.
     */
    void ReMesh();
    
    /**
     * Check consistency of our internal data structures.
     */
    void Validate();
    
    std::vector<Node<DIM> >& rGetNodes();
    const std::vector<Node<DIM> >& rGetNodes() const;
    
    /**
     * Move a cell to a new location.
     * @param iter  pointer to the cell to move
     * @param rNewLocation  where to move it to
     */
    void MoveCell(typename AbstractTissue<DIM>::Iterator iter, ChastePoint<DIM>& rNewLocation);
    
};

template<unsigned DIM>
SimpleTissue<DIM>::SimpleTissue(const std::vector<Node<DIM> >& rNodes, 
                                const std::vector<TissueCell>& rCells)
        : AbstractTissue<DIM>(rCells),
          mNodes(rNodes)
{
    Validate();    
}

template<unsigned DIM>
SimpleTissue<DIM>::SimpleTissue(const std::vector<Node<DIM> >& rNodes)
        : AbstractTissue<DIM>(),
          mNodes(rNodes)
{
}


template<unsigned DIM>
void SimpleTissue<DIM>::Validate()
{
    std::vector<bool> validated_node(GetNumNodes()); 
    
    for (typename AbstractTissue<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        unsigned node_index = cell_iter->GetNodeIndex();
        validated_node[node_index] = true;
    }
    
    for (unsigned i=0; i<validated_node.size(); i++)
    {
        if(!validated_node[i])
        {
            std::stringstream ss;
            ss << "Node " << i << " does not appear to have a cell associated with it";
            EXCEPTION(ss.str()); 
        }
    }
}

template<unsigned DIM>
std::vector<Node<DIM> >& SimpleTissue<DIM>::rGetNodes()
{
    return mNodes;
}

template<unsigned DIM>
const std::vector<Node<DIM> >& SimpleTissue<DIM>::rGetNodes() const
{
    return mNodes;
}

template<unsigned DIM>
Node<DIM>* SimpleTissue<DIM>::GetNode(unsigned index)
{
    return &(mNodes[index]);
}

template<unsigned DIM>
void SimpleTissue<DIM>::SetNode(unsigned index, ChastePoint<DIM> point)
{
    mNodes[index].SetPoint(point);    
}

template<unsigned DIM>
void SimpleTissue<DIM>::MoveCell(typename AbstractTissue<DIM>::Iterator iter, ChastePoint<DIM>& rNewLocation)
{
    unsigned index = iter.GetNode()->GetIndex();
    SetNode(index, rNewLocation);
}

template<unsigned DIM>
TissueCell* SimpleTissue<DIM>::AddCell(TissueCell newCell, c_vector<double,DIM> newLocation)
{
    // Create a new node
    Node<DIM> new_node(GetNumNodes(), newLocation, false); // never on boundary
    
    unsigned new_node_index = AddNode(&new_node); //Uses copy constructor (so it doesn't matter that new_node goes out of scope)

    // Associate the new cell with the node
    newCell.SetNodeIndex(new_node_index);
    this->mCells.push_back(newCell);
    
    TissueCell *p_created_cell = &(this->mCells.back());
    this->mNodeCellMap[new_node_index] = p_created_cell;
    
    return p_created_cell;
}

template<unsigned DIM>
void SimpleTissue<DIM>::ReMesh()
{
    // Create and reserve space for a temporary vector
    /// \todo: reserve space equal to mNodes.size() for this vector (see #642)
    std::vector<Node<DIM> > old_nodes;
    
    // Store all non-deleted nodes in the temporary vector
    for (unsigned i=0; i<mNodes.size(); i++)
    {
        if (mNodes[i].IsDeleted()==false)
        {
            old_nodes.push_back(mNodes[i]);
        }
    }
    
    // Update mNodes
    mNodes = old_nodes;
    
    // Update the correspondence between nodes and cells.    
    // We expect the node indices to be {0,1,...,num_nodes}
    std::set<unsigned> expected_node_indices;    
    for (unsigned i=0; i<GetNumNodes(); i++)
    {
        expected_node_indices.insert(i);
    }
    
    // Get the actual set of node indices
    std::set<unsigned> node_indices;
    for (typename AbstractTissue<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        unsigned node_index = cell_iter->GetNodeIndex();
        node_indices.insert(node_index);       
    }
    
    // If necessary, update the node cell map
    if (node_indices != expected_node_indices)
    {
        // Fix up the mappings between cells and nodes
        this->mNodeCellMap.clear();
        unsigned new_node_index = 0;
        for (typename AbstractTissue<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
        {
            cell_iter->SetNodeIndex(new_node_index);
            this->mNodeCellMap[new_node_index] = &(*cell_iter);
            new_node_index++;
        }
    }
    
    Validate();
}

template<unsigned DIM>
unsigned SimpleTissue<DIM>::RemoveDeadCells()
{
    unsigned num_removed = 0;
    
    for (std::list<TissueCell>::iterator cell_iter = this->mCells.begin();
         cell_iter != this->mCells.end();
         ++cell_iter)
    {
        if (cell_iter->IsDead())
        {   
            // Remove the node from the mesh
            num_removed++;
            this->GetNodeCorrespondingToCell(*cell_iter)->MarkAsDeleted();
            cell_iter = this->mCells.erase(cell_iter);
            --cell_iter;
        }
    }
    return num_removed;
}

template<unsigned DIM>
unsigned SimpleTissue<DIM>::AddNode(Node<DIM> *pNewNode)
{
    /// \todo: employ a std::vector of deleted node indices to re-use indices? 
    pNewNode->SetIndex(mNodes.size());
    mNodes.push_back(*pNewNode);        
    return pNewNode->GetIndex();
}    

template<unsigned DIM>
unsigned SimpleTissue<DIM>::GetNumNodes()
{
    return mNodes.size();
}


#include "TemplatedExport.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SimpleTissue)

namespace boost {
namespace serialization {

/**
 * Non-intrusive serialization for Node - save method.
 */
template<class Archive, unsigned SPACE_DIM>
inline void save(
    Archive & ar,
    const Node<SPACE_DIM> &rNode,
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
    Node<SPACE_DIM> &rNode,
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
    
    // Invoke inplace constructor to initialize instance
    ::new(t)Node<DIM>(index, loc, is_boundary);
}


/**
 * Serialize information required to construct a SimpleTissue facade.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const SimpleTissue<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    ar & t->rGetNodes();
}

/**
 * De-serialize constructor parameters and initialise SimpleTissue.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, SimpleTissue<DIM> * t, const unsigned int file_version)
{
    // Load the nodes
    std::vector<Node<DIM> > nodes;
    ar >> nodes;
    
    // Invoke inplace constructor to initialize instance
    ::new(t)SimpleTissue<DIM>(nodes);
}

}} // close namespaces


#endif /*SIMPLETISSUE_HPP_*/
