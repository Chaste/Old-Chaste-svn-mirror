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

#include <map>
#include "NodesOnlyMesh.hpp"

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::ConstructNodesWithoutMesh(const std::vector<Node<SPACE_DIM>*>& rNodes)
{
    this->Clear();
    for (unsigned i=0; i<rNodes.size(); i++)
    {
        assert(!rNodes[i]->IsDeleted());
        bool boundary = rNodes[i]->IsBoundaryNode();
        c_vector<double, SPACE_DIM> location = rNodes[i]->rGetLocation();
        
        Node<SPACE_DIM>* p_node_copy = new Node<SPACE_DIM>(i, location, boundary);
        this->mNodes.push_back(p_node_copy);
        mCellRadii.push_back(1.0);
    }
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::ConstructNodesWithoutMesh(const TetrahedralMesh<SPACE_DIM,SPACE_DIM>& rGeneratingMesh)
{
    ConstructNodesWithoutMesh(rGeneratingMesh.mNodes);
}

template<unsigned SPACE_DIM>
double NodesOnlyMesh<SPACE_DIM>::GetCellRadius(unsigned index)
{
    return mCellRadii[index];    
}
    
template<unsigned SPACE_DIM>     
void NodesOnlyMesh<SPACE_DIM>::SetCellRadius(unsigned index, double radius)
{
    mCellRadii[index] = radius;
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::ReMesh(NodeMap& map)
{
    // Store the node locations
    std::vector<c_vector<double, SPACE_DIM> > old_node_locations;
    unsigned new_index = 0;
    for (unsigned i=0; i<this->GetNumAllNodes(); i++)
    {
        if (this->mNodes[i]->IsDeleted())
        {
            map.SetDeleted(i);
        }
        else
        {
            map.SetNewIndex(i, new_index);
            old_node_locations.push_back(this->mNodes[i]->rGetLocation());
            new_index++;
        }
    }

    // Remove current data
    this->Clear();

    // Construct the nodes and boundary nodes
    for (unsigned node_index=0; node_index<old_node_locations.size(); node_index++)
    {
        Node<SPACE_DIM>* p_node = new Node<SPACE_DIM>(node_index, old_node_locations[node_index], false);
        this->mNodes.push_back(p_node);

        // As we're in 1D, the boundary nodes are simply at either end of the mesh
        if (node_index==0 || node_index==old_node_locations.size()-1)
        {
            this->mBoundaryNodes.push_back(p_node);
        }
    }

    // Create a map between node indices and node locations
    std::map<double, unsigned> location_index_map;
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        location_index_map[this->mNodes[i]->rGetLocation()[0]] = this->mNodes[i]->GetIndex();
    }

    // Use this map to generate a vector of node indices that are ordered spatially
    std::vector<unsigned> node_indices_ordered_spatially;
    for (std::map<double, unsigned>::iterator iter = location_index_map.begin();
         iter != location_index_map.end();
         ++iter)
    {
        node_indices_ordered_spatially.push_back(iter->second);
    }
}

template<unsigned SPACE_DIM>
void NodesOnlyMesh<SPACE_DIM>::DeleteNode(unsigned index)
{
    if (this->mNodes[index]->IsDeleted())
    {
        EXCEPTION("Trying to delete a deleted node");
    }
    
    this->mNodes[index]->MarkAsDeleted();
    this->mDeletedNodeIndices.push_back(index);
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class NodesOnlyMesh<1>;
template class NodesOnlyMesh<2>;
template class NodesOnlyMesh<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NodesOnlyMesh)
