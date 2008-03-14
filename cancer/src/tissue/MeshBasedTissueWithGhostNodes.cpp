#ifndef MESHBASEDTISSUEWITHGHOSTNODES_CPP
#define MESHBASEDTISSUEWITHGHOSTNODES_CPP

#include "MeshBasedTissueWithGhostNodes.hpp"

///\todo: Make this constructor take in ghost nodes, and validate the three objects
// are in sync ie num cells + num ghost nodes = num_nodes ? this would mean all ghosts
// *cannot* be cells, making it more difficult to construct the cells.
// also check cell.GetNodeIndices() is in the mesh, and covers the mesh, etc.
template<unsigned DIM>
MeshBasedTissueWithGhostNodes<DIM>::MeshBasedTissueWithGhostNodes(
     ConformingTetrahedralMesh<DIM, DIM>& rMesh,
     const std::vector<TissueCell>& rCells,
     const std::set<unsigned> ghostNodeIndices,
     bool deleteMesh)
             : MeshBasedTissue<DIM>(rMesh, rCells, deleteMesh)
{
    SetGhostNodes(ghostNodeIndices);
    this->mTissueContainsGhostNodes = true;
    ValidateWithGhostNodes();
}

template<unsigned DIM>
MeshBasedTissueWithGhostNodes<DIM>::MeshBasedTissueWithGhostNodes(ConformingTetrahedralMesh<DIM, DIM>& rMesh)
             : MeshBasedTissue<DIM>(rMesh)
{
    this->mTissueContainsGhostNodes = true;
}

template<unsigned DIM>
std::vector<bool>& MeshBasedTissueWithGhostNodes<DIM>::rGetGhostNodes()
{
    return this->mIsGhostNode;
}

template<unsigned DIM>
unsigned MeshBasedTissueWithGhostNodes<DIM>::GetGhostNodesSize()
{
    return this->mIsGhostNode.size();
}

template<unsigned DIM>
bool MeshBasedTissueWithGhostNodes<DIM>::IsGhostNode(unsigned index)
{
    return this->mIsGhostNode[index];
}

template<unsigned DIM>
std::set<unsigned> MeshBasedTissueWithGhostNodes<DIM>::GetGhostNodeIndices()
{
    std::set<unsigned> ghost_node_indices;
    for (unsigned i=0; i<this->mIsGhostNode.size(); i++)
    {
        if (this->mIsGhostNode[i])
        {
            ghost_node_indices.insert(i);    
        }        
    }
    return ghost_node_indices;        
}

template<unsigned DIM>
void MeshBasedTissueWithGhostNodes<DIM>::SetGhostNodes(const std::vector<bool>& rGhostNodes)
{
    this->mIsGhostNode = rGhostNodes;
}

template<unsigned DIM> 
void MeshBasedTissueWithGhostNodes<DIM>::SetGhostNodes(const std::set<unsigned>& ghostNodeIndices)
{
    // Reinitialise all entries of mIsGhostNode to false
    this->mIsGhostNode = std::vector<bool>(this->mrMesh.GetNumNodes(), false);

    // Update mIsGhostNode
    std::set<unsigned>::iterator iter = ghostNodeIndices.begin();
    while(iter!=ghostNodeIndices.end())
    {
        this->mIsGhostNode[*iter] = true;
        iter++;
    }
}

template<unsigned DIM>
void MeshBasedTissueWithGhostNodes<DIM>::UpdateGhostPositions(double dt)
{
    // Initialise vector of ghost node velocities
    std::vector<c_vector<double, DIM> > drdt(this->GetNumNodes());
    for (unsigned i=0; i<drdt.size(); i++)
    {
        drdt[i] = zero_vector<double>(DIM);
    }

    // Calculate ghost node velocities
    for (typename ConformingTetrahedralMesh<DIM, DIM>::EdgeIterator edge_iterator=this->mrMesh.EdgesBegin();
        edge_iterator!=this->mrMesh.EdgesEnd();
        ++edge_iterator)
    {
        unsigned nodeA_global_index = edge_iterator.GetNodeA()->GetIndex();
        unsigned nodeB_global_index = edge_iterator.GetNodeB()->GetIndex();

        c_vector<double, DIM> force = CalculateForceBetweenNodes(nodeA_global_index, nodeB_global_index);
         
        double damping_constant = CancerParameters::Instance()->GetDampingConstantNormal();
                
        if (!this->mIsGhostNode[nodeA_global_index])
        {
            drdt[nodeB_global_index] -= force / damping_constant;
        }
        else
        {
            drdt[nodeA_global_index] += force / damping_constant;
                
            if (this->mIsGhostNode[nodeB_global_index])
            {
                drdt[nodeB_global_index] -= force / damping_constant;
            }
        }
    }
    
    for (unsigned index = 0; index<this->GetNumNodes(); index++)
    {
        if ((!this->GetNode(index)->IsDeleted()) && this->mIsGhostNode[index])
        {
            ChastePoint<DIM> new_point(this->GetNode(index)->rGetLocation() + dt*drdt[index]);
            this->mrMesh.SetNode(index, new_point, false);
        }
    }
}

/**
 * Calculates the force between two nodes.
 * 
 * Note that this assumes they are connected
 * 
 * @param NodeAGlobalIndex
 * @param NodeBGlobalIndex
 * 
 * @return The force exerted on Node A by Node B.
 */
template<unsigned DIM> 
c_vector<double, DIM> MeshBasedTissueWithGhostNodes<DIM>::CalculateForceBetweenNodes(const unsigned& rNodeAGlobalIndex, const unsigned& rNodeBGlobalIndex)
{
    assert(rNodeAGlobalIndex!=rNodeBGlobalIndex);
    c_vector<double, DIM> unit_difference;
    c_vector<double, DIM> node_a_location = this->mrMesh.GetNode(rNodeAGlobalIndex)->rGetLocation();
    c_vector<double, DIM> node_b_location = this->mrMesh.GetNode(rNodeBGlobalIndex)->rGetLocation();
    
    // There is reason not to substract one position from the other (cylindrical meshes)
    unit_difference = this->mrMesh.GetVectorFromAtoB(node_a_location, node_b_location);   
    
    double distance_between_nodes = norm_2(unit_difference);
    
    unit_difference /= distance_between_nodes;
    
    double rest_length = 1.0;
    
    return CancerParameters::Instance()->GetSpringStiffness() * unit_difference * (distance_between_nodes - rest_length);
}

template<unsigned DIM>  
TissueCell* MeshBasedTissueWithGhostNodes<DIM>::AddCell(TissueCell newCell, c_vector<double,DIM> newLocation)
{
    // Add cell
    TissueCell *p_created_cell = MeshBasedTissue<DIM>::AddCell(newCell, newLocation);

    // Update size of mIsGhostNode if necessary    
    unsigned new_node_index = p_created_cell->GetNodeIndex();
    
    if (this->GetNumNodes() > this->mIsGhostNode.size())
    {
        this->mIsGhostNode.resize(this->GetNumNodes());
        this->mIsGhostNode[new_node_index] = false;
    }
    return p_created_cell;
}

template<unsigned DIM>
void MeshBasedTissueWithGhostNodes<DIM>::ValidateWithGhostNodes()
{    
    std::vector<bool> validated_node = mIsGhostNode; 
    for (typename AbstractTissue<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        unsigned node_index = cell_iter->GetNodeIndex();
        validated_node[node_index] = true;
    }
    
    for (unsigned i=0; i<validated_node.size(); i++)
    {
        if (!validated_node[i])
        {
            std::stringstream ss;
            ss << "Node " << i << " does not appear to be a ghost node or have a cell associated with it";
            EXCEPTION(ss.str());
        }
    }
}

template<unsigned DIM>
void MeshBasedTissueWithGhostNodes<DIM>::UpdateGhostNodesAfterReMesh(NodeMap& rMap)
{
    // Copy mIsGhostNode to a temporary vector
    std::vector<bool> ghost_nodes_before_remesh = mIsGhostNode;
    
    // Reinitialise mIsGhostNode
    mIsGhostNode.clear();
    mIsGhostNode.resize(this->GetNumNodes());

    // Update mIsGhostNode using the node map
    for (unsigned old_index=0; old_index<rMap.Size(); old_index++)
    {
        if (!rMap.IsDeleted(old_index))
        {
            unsigned new_index = rMap.GetNewIndex(old_index);
            mIsGhostNode[new_index] = ghost_nodes_before_remesh[old_index];
        }
    }
}

#endif //MESHBASEDTISSUEWITHGHOSTNODES_CPP

