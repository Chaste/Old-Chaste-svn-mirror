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
#include "Debug.hpp"

template<unsigned DIM>
NodeBasedTissue<DIM>::NodeBasedTissue(const std::vector<Node<DIM>* > nodes,
                                      const std::vector<TissueCell>& rCells,
                                      const std::vector<unsigned> locationIndices,
                                      bool deleteNodes)
    : AbstractCellCentreBasedTissue<DIM>(rCells, locationIndices),
      mNodes(nodes.begin(), nodes.end()),
      mDeleteNodes(deleteNodes)
{
    Clear();
    mAddedNodes = true;
    Validate();
    
}


template<unsigned DIM>
NodeBasedTissue<DIM>::NodeBasedTissue(const std::vector<Node<DIM>* > nodes, bool deleteNodes)
    : AbstractCellCentreBasedTissue<DIM>(),
      mNodes(nodes.begin(), nodes.end()),
      mDeleteNodes(deleteNodes)
{
    Clear();
    mAddedNodes = true;
}


template<unsigned DIM>
NodeBasedTissue<DIM>::NodeBasedTissue(const AbstractMesh<DIM,DIM>& rMesh,
                                      const std::vector<TissueCell>& rCells)
    : AbstractCellCentreBasedTissue<DIM>(rCells),
      mDeleteNodes(true)
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
    if (mDeleteNodes)
    {
        for (unsigned i=0; i<mNodes.size(); i++)
        {
            delete mNodes[i];
        }
    }
}


template<unsigned DIM>
void NodeBasedTissue<DIM>::Clear()
{
    mpNodeBoxCollection = NULL;
    mNodePairs.clear();
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
void NodeBasedTissue<DIM>::SplitUpIntoBoxes(double cutOffLength, c_vector<double, 2*DIM> domainSize)
{
    mpNodeBoxCollection = new NodeBoxCollection<DIM>(cutOffLength, domainSize);

    for (unsigned i=0; i<mNodes.size(); i++)
    {
        unsigned box_index = mpNodeBoxCollection->CalculateContainingBox(mNodes[i]);
        mpNodeBoxCollection->rGetBox(box_index).AddNode(mNodes[i]);
    }       
}

template<unsigned DIM>
void NodeBasedTissue<DIM>::FindMaxAndMin()
{
    c_vector<double, DIM> min_posn;
    c_vector<double, DIM> max_posn;
    for (unsigned i=0; i<DIM; i++)
    {
        min_posn(i) = DBL_MAX;
        max_posn(i) = -DBL_MIN;
    }
    for (unsigned i=0; i<mNodes.size(); i++)
    {
        for (unsigned j=0; j<DIM; j++)
        {
            if (this->GetNode(i)->rGetLocation()[j] > max_posn(j))
            {
                max_posn(j) = this->GetNode(i)->rGetLocation()[j];
            }
            if (this->GetNode(i)->rGetLocation()[j] < min_posn(j))
            {
                min_posn(j) = this->GetNode(i)->rGetLocation()[j];
            }
        }
    }
    
    for (unsigned i=0; i<DIM; i++)
    {
        assert(min_posn(i) != DBL_MAX);
        mMinSpatialPositions(i) = min_posn(i);
        
        assert(max_posn(i) != -DBL_MIN);
        mMaxSpatialPositions(i) = max_posn(i);
    }
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
void NodeBasedTissue<DIM>::Update(bool hasHadBirthsOrDeaths)
{
    if (hasHadBirthsOrDeaths)
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
    if (mpNodeBoxCollection==NULL)
    {
        FindMaxAndMin();
        // Something here to set up the domain size (max and min of each node position dimension)
        c_vector<double, 2*DIM> domainSize;
        
        for (unsigned i=0; i<DIM; i++)
        {
            domainSize(2*i) = mMinSpatialPositions(i);
            domainSize(2*i+1) = mMaxSpatialPositions(i);
        }
        
        // Add this cancer parameter and suggest that mechanics systems set it.
        SplitUpIntoBoxes(CancerParameters::Instance()->GetMechanicsCutOffLength(), domainSize);
    }
    mpNodeBoxCollection->CalculateNodePairs(mNodes,mNodePairs);
    /// \todo #899 update the node pairs here.
//    c_vector<double, 2*DIM>& mpNodeBoxCollection->rGetBox(CalculateContainingBox(this->GetNode(node_index)))
//                                                                                            ->rGetMinAndMaxValues();
//        
//        for (unsigned i=0; i<DIM; i++)
//        {
//            if (new_node_location[i] < )
//        c_vector<double, 2*DIM>& rGetMinAndMaxValues();
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



template<unsigned DIM>
NodeBoxCollection<DIM>* NodeBasedTissue<DIM>::GetNodeBoxCollection()
{
    return mpNodeBoxCollection;   
}

template<unsigned DIM>
std::set< std::pair<Node<DIM>*, Node<DIM>* > > NodeBasedTissue<DIM>::GetNodePairs()
{
    return mNodePairs;
}
/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////


template class NodeBasedTissue<1>;
template class NodeBasedTissue<2>;
template class NodeBasedTissue<3>;
