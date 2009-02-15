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
#include "NodeBasedTissue.hpp"

template<unsigned DIM>
NodeBasedTissue<DIM>::NodeBasedTissue(const std::vector<Node<DIM>* > nodes,
                                      const std::vector<TissueCell>& rCells,
                                      const std::vector<unsigned> locationIndices)
    : AbstractCellCentreBasedTissue<DIM>(rCells, locationIndices),
      mNodes(nodes.begin(), nodes.end())
{
    Clear();
    mAddedNodes = true;
    Validate();
}


template<unsigned DIM>
NodeBasedTissue<DIM>::NodeBasedTissue(const std::vector<Node<DIM>* > nodes)
    : AbstractCellCentreBasedTissue<DIM>(),
      mNodes(nodes.begin(), nodes.end())
{
    Clear();
    mAddedNodes = true;
}


template<unsigned DIM>
NodeBasedTissue<DIM>::NodeBasedTissue(const AbstractMesh<DIM,DIM>& rMesh,
                                      const std::vector<TissueCell>& rCells)
    : AbstractCellCentreBasedTissue<DIM>(rCells)
{
    Clear();
    mNodes.reserve(rMesh.GetNumNodes());
    // Copy the actual node objects
    for (unsigned i=0; i<rMesh.GetNumNodes(); i++)
    {
        Node<DIM>* p_node = new Node<DIM>(*(rMesh.GetNode(i)));
        mNodes.push_back(p_node);
    }
    mAddedNodes = true;
    Validate();
}


template<unsigned DIM>
NodeBasedTissue<DIM>::~NodeBasedTissue()
{
    Clear();
    // Free node memory
    for (unsigned i=0; i<mNodes.size(); i++)
    {
        delete mNodes[i];
    }
}


template<unsigned DIM>
void NodeBasedTissue<DIM>::Clear()
{
    mDeletedNodeIndices.clear();
    mAddedNodes = false;
}


template<unsigned DIM>
void NodeBasedTissue<DIM>::Validate()
{
    std::vector<bool> validated_node(GetNumNodes());

    for (typename AbstractTissue<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        unsigned node_index = this->mCellLocationMap[&(*cell_iter)];
        validated_node[node_index] = true;
    }

    for (unsigned i=0; i<validated_node.size(); i++)
    {
        if (!validated_node[i])
        {
            std::stringstream ss;
            ss << "Node " << i << " does not appear to have a cell associated with it";
            EXCEPTION(ss.str());
        }
    }
}


template<unsigned DIM>
std::vector<Node<DIM>* >& NodeBasedTissue<DIM>::rGetNodes()
{
    return mNodes;
}


template<unsigned DIM>
const std::vector<Node<DIM>* >& NodeBasedTissue<DIM>::rGetNodes() const
{
    return mNodes;
}


template<unsigned DIM>
Node<DIM>* NodeBasedTissue<DIM>::GetNode(unsigned index)
{
    return mNodes[index];
}


template<unsigned DIM>
void NodeBasedTissue<DIM>::SetNode(unsigned nodeIndex, ChastePoint<DIM>& rNewLocation)
{
    mNodes[nodeIndex]->SetPoint(rNewLocation);
}


template<unsigned DIM>
void NodeBasedTissue<DIM>::Update()
{
    // Create and reserve space for a temporary vector
    std::vector<Node<DIM>* > old_nodes;
    old_nodes.reserve(mNodes.size());

    // Store all non-deleted nodes in the temporary vector
    for (unsigned i=0; i<mNodes.size(); i++)
    {
        if ( !mNodes[i]->IsDeleted() )
        {
            old_nodes.push_back(mNodes[i]);
        }
        else
        {
            // Free node memory
            delete mNodes[i];
        }
    }

    std::map<unsigned,TissueCell*> old_map = this->mLocationCellMap;
    mNodes.clear();
    
    // Clear maps
    this->mLocationCellMap.clear();
    this->mCellLocationMap.clear();

    // Update mNodes to new indices which go from 0 to NumNodes-1.
    for (unsigned i=0; i<old_nodes.size() ; i++)
    {
        // Get the living cell associated with the old node
        TissueCell* p_live_cell = old_map[old_nodes[i]->GetIndex()];
        // Set the node up
        mNodes.push_back(old_nodes[i]);
        mNodes[i]->SetIndex(i);
        // Set the maps up
        this->mLocationCellMap[i] = p_live_cell;
        this->mCellLocationMap[p_live_cell] = i;
    }


    // Remove current dead indices data
    Clear();

    Validate();
}


template<unsigned DIM>
unsigned NodeBasedTissue<DIM>::RemoveDeadCells()
{
    unsigned num_removed = 0;

    for (std::list<TissueCell>::iterator cell_iter = this->mCells.begin();
         cell_iter != this->mCells.end();
         ++cell_iter)
    {
        if (cell_iter->IsDead())
        {
            // Remove the node
            num_removed++;
            this->GetNodeCorrespondingToCell(&(*cell_iter))->MarkAsDeleted();
            mDeletedNodeIndices.push_back( this->mCellLocationMap[&(*cell_iter)] );
            cell_iter = this->mCells.erase(cell_iter);
            --cell_iter;
        }
    }
    return num_removed;
}


template<unsigned DIM>
unsigned NodeBasedTissue<DIM>::AddNode(Node<DIM>* pNewNode)
{
    if (mDeletedNodeIndices.empty())
    {
        pNewNode->SetIndex(mNodes.size());
        mNodes.push_back(pNewNode);
    }
    else
    {
        unsigned index = mDeletedNodeIndices.back();
        pNewNode->SetIndex(index);
        mDeletedNodeIndices.pop_back();
        delete mNodes[index];
        mNodes[index] = pNewNode;
    }
    mAddedNodes = true;
    return pNewNode->GetIndex();
}


template<unsigned DIM>
unsigned NodeBasedTissue<DIM>::GetNumNodes()
{
    return mNodes.size() - mDeletedNodeIndices.size();
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////


template class NodeBasedTissue<1>;
template class NodeBasedTissue<2>;
template class NodeBasedTissue<3>;
