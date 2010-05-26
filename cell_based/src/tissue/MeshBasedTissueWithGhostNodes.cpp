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
#include "MeshBasedTissueWithGhostNodes.hpp"


template<unsigned DIM>
MeshBasedTissueWithGhostNodes<DIM>::MeshBasedTissueWithGhostNodes(
     MutableMesh<DIM, DIM>& rMesh,
     std::vector<TissueCell>& rCells,
     const std::vector<unsigned> locationIndices,
     bool deleteMesh)
             : MeshBasedTissue<DIM>(rMesh, rCells, locationIndices, deleteMesh, false) // do not call the base class Validate()
{
    if (!locationIndices.empty())
    {
        // Create a set of node indices corresponding to ghost nodes
        std::set<unsigned> node_indices;
        std::set<unsigned> location_indices;
        std::set<unsigned> ghost_node_indices;

        for (unsigned i=0; i<this->GetNumNodes(); i++)
        {
            node_indices.insert(this->GetNode(i)->GetIndex());
        }
        for (unsigned i=0; i<locationIndices.size(); i++)
        {
            location_indices.insert(locationIndices[i]);
        }

        std::set_difference(node_indices.begin(), node_indices.end(),
                            location_indices.begin(), location_indices.end(),
                            std::inserter(ghost_node_indices, ghost_node_indices.begin()));

        // This method finishes and then calls Validate()
        SetGhostNodes(ghost_node_indices);
    }
    else
    {
        this->mIsGhostNode = std::vector<bool>(this->GetNumNodes(), false);
        Validate();
    }
}

template<unsigned DIM>
MeshBasedTissueWithGhostNodes<DIM>::MeshBasedTissueWithGhostNodes(MutableMesh<DIM, DIM>& rMesh)
             : MeshBasedTissue<DIM>(rMesh)
{
}

template<unsigned DIM>
std::vector<bool>& MeshBasedTissueWithGhostNodes<DIM>::rGetGhostNodes()
{
    return this->mIsGhostNode;
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
void MeshBasedTissueWithGhostNodes<DIM>::SetGhostNodes(const std::set<unsigned>& rGhostNodeIndices)
{
    // Reinitialise all entries of mIsGhostNode to false
    this->mIsGhostNode = std::vector<bool>(this->mrMesh.GetNumNodes(), false);

    // Update mIsGhostNode
    for (std::set<unsigned>::iterator iter=rGhostNodeIndices.begin(); iter!=rGhostNodeIndices.end(); ++iter)
    {
        this->mIsGhostNode[*iter] = true;
    }

    Validate();
}

template<unsigned DIM>
void MeshBasedTissueWithGhostNodes<DIM>::UpdateGhostPositions(double dt)
{
    // Initialise vector of forces on ghost nodes
    std::vector<c_vector<double, DIM> > drdt(this->GetNumNodes());
    for (unsigned i=0; i<drdt.size(); i++)
    {
        drdt[i] = zero_vector<double>(DIM);
    }

    // Calculate forces on ghost nodes
    for (typename MutableMesh<DIM, DIM>::EdgeIterator edge_iterator=this->mrMesh.EdgesBegin();
        edge_iterator!=this->mrMesh.EdgesEnd();
        ++edge_iterator)
    {
        unsigned nodeA_global_index = edge_iterator.GetNodeA()->GetIndex();
        unsigned nodeB_global_index = edge_iterator.GetNodeB()->GetIndex();

        c_vector<double, DIM> force = CalculateForceBetweenNodes(nodeA_global_index, nodeB_global_index);

        double damping_constant = TissueConfig::Instance()->GetDampingConstantNormal();

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

    for (typename AbstractMesh<DIM,DIM>::NodeIterator node_iter = this->mrMesh.GetNodeIteratorBegin();
         node_iter != this->mrMesh.GetNodeIteratorEnd();
         ++node_iter)
    {
        unsigned node_index = node_iter->GetIndex();
        if (this->mIsGhostNode[node_index])
        {
            ChastePoint<DIM> new_point(node_iter->rGetLocation() + dt*drdt[node_index]);
            this->mrMesh.SetNode(node_index, new_point, false);
        }
    }
}

template<unsigned DIM>
c_vector<double, DIM> MeshBasedTissueWithGhostNodes<DIM>::CalculateForceBetweenNodes(const unsigned& rNodeAGlobalIndex, const unsigned& rNodeBGlobalIndex)
{
    assert(rNodeAGlobalIndex!=rNodeBGlobalIndex);
    c_vector<double, DIM> unit_difference;
    c_vector<double, DIM> node_a_location = this->GetNode(rNodeAGlobalIndex)->rGetLocation();
    c_vector<double, DIM> node_b_location = this->GetNode(rNodeBGlobalIndex)->rGetLocation();

    // There is reason not to subtract one position from the other (cylindrical meshes)
    unit_difference = this->mrMesh.GetVectorFromAtoB(node_a_location, node_b_location);

    double distance_between_nodes = norm_2(unit_difference);

    unit_difference /= distance_between_nodes;

    double rest_length = 1.0;

    return TissueConfig::Instance()->GetMeinekeSpringStiffness() * unit_difference * (distance_between_nodes - rest_length);
}

template<unsigned DIM>
TissueCell* MeshBasedTissueWithGhostNodes<DIM>::AddCell(TissueCell& rNewCell, const c_vector<double,DIM>& rCellDivisionVector, TissueCell* pParentCell)
{
    // Add new cell to tissue
    TissueCell* p_created_cell = MeshBasedTissue<DIM>::AddCell(rNewCell, rCellDivisionVector, pParentCell);

    // Update size of mIsGhostNode if necessary
    unsigned new_node_index = this->mCellLocationMap[p_created_cell];

    if (this->GetNumNodes() > this->mIsGhostNode.size())
    {
        this->mIsGhostNode.resize(this->GetNumNodes());
        this->mIsGhostNode[new_node_index] = false;
    }

    // Return pointer to new cell
    return p_created_cell;
}

template<unsigned DIM>
void MeshBasedTissueWithGhostNodes<DIM>::Validate()
{
    // Get a list of all the nodes that are ghosts
    std::vector<bool> validated_node = mIsGhostNode;
    assert(mIsGhostNode.size()==this->GetNumNodes());

    // Look through all of the cells and record what node they are associated with.
    for (typename AbstractTissue<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        unsigned node_index = this->mCellLocationMap[&(*cell_iter)];

        // If the node attached to this cell is labelled as a ghost node, then throw an error
        if (mIsGhostNode[node_index])
        {
            std::stringstream ss;
            ss << "Node " << node_index << " is labelled as a ghost node and has a cell attached";
            EXCEPTION(ss.str());
        }
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

template<unsigned DIM>
void MeshBasedTissueWithGhostNodes<DIM>::UpdateNodeLocations(const std::vector< c_vector<double, DIM> >& rNodeForces, double dt)
{
    // First update ghost positions first because they do not affect the real cells
    UpdateGhostPositions(dt);

    // Then call the base class method
    AbstractCellCentreBasedTissue<DIM>::UpdateNodeLocations(rNodeForces, dt);
}

template<unsigned DIM>
void MeshBasedTissueWithGhostNodes<DIM>::WriteVtkResultsToFile()
{
#ifdef CHASTE_VTK
    VertexMeshWriter<DIM, DIM> mesh_writer(this->mDirPath, "results", false);
    std::stringstream time;
    time << SimulationTime::Instance()->GetTimeStepsElapsed();

    unsigned num_elements = this->mpVoronoiTessellation->GetNumElements();
    std::vector<double> ghosts(num_elements);
    std::vector<double> cell_types(num_elements);
    std::vector<double> cell_ancestors(num_elements);
    std::vector<double> cell_mutation_states(num_elements);
    std::vector<double> cell_ages(num_elements);
    std::vector<double> cell_cycle_phases(num_elements);
    std::vector<double> cell_volumes(num_elements);

    // Loop over Voronoi elements
    for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = this->mpVoronoiTessellation->GetElementIteratorBegin();
         elem_iter != this->mpVoronoiTessellation->GetElementIteratorEnd();
         ++elem_iter)
    {
        // Get index of this element in the Voronoi tessellation mesh
        unsigned elem_index = elem_iter->GetIndex();

        unsigned node_index = this->mpVoronoiTessellation->GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex(elem_index);

        ghosts[elem_index] = (double)(this->IsGhostNode(node_index));

        if (!this->IsGhostNode(node_index))
        {
            // Get the cell corresponding to this element
            TissueCell* p_cell = this->mLocationCellMap[node_index];

            if (TissueConfig::Instance()->GetOutputCellAncestors())
            {
                double ancestor_index = (p_cell->GetAncestor() == UNSIGNED_UNSET) ? (-1.0) : (double)p_cell->GetAncestor();
                cell_ancestors[elem_index] = ancestor_index;
            }
            if (TissueConfig::Instance()->GetOutputCellProliferativeTypes())
            {
                double cell_type = p_cell->GetCellCycleModel()->GetCellProliferativeType();
                cell_types[elem_index] = cell_type;
            }
            if (TissueConfig::Instance()->GetOutputCellMutationStates())
            {
                double mutation_state = p_cell->GetMutationState()->GetColour();
                cell_mutation_states[elem_index] = mutation_state;
            }
            if (TissueConfig::Instance()->GetOutputCellAges())
            {
                double age = p_cell->GetAge();
                cell_ages[elem_index] = age;
            }
            if (TissueConfig::Instance()->GetOutputCellCyclePhases())
            {
                double cycle_phase = p_cell->GetCellCycleModel()->GetCurrentCellCyclePhase();
                cell_cycle_phases[elem_index] = cycle_phase;
            }
            if (TissueConfig::Instance()->GetOutputCellVolumes())
            {
                double cell_volume = this->mpVoronoiTessellation->GetVolumeOfElement(elem_index);
                cell_volumes[elem_index] = cell_volume;
            }
        }
        else
        {
            if (TissueConfig::Instance()->GetOutputCellAncestors())
            {
                cell_ancestors[elem_index] = -1.0;
            }
            if (TissueConfig::Instance()->GetOutputCellProliferativeTypes())
            {
                cell_types[elem_index] = -1.0;
            }
            if (TissueConfig::Instance()->GetOutputCellMutationStates())
            {
                cell_mutation_states[elem_index] = -1.0;
            }
            if (TissueConfig::Instance()->GetOutputCellAges())
            {
                cell_ages[elem_index] = -1.0;
            }
            if (TissueConfig::Instance()->GetOutputCellCyclePhases())
            {
                cell_cycle_phases[elem_index] = -1.0;
            }
            if (TissueConfig::Instance()->GetOutputCellVolumes())
            {
                cell_volumes[elem_index] = -1.0;
            }
        }
    }

    mesh_writer.AddCellData("Non-ghosts", ghosts);
    if (TissueConfig::Instance()->GetOutputCellProliferativeTypes())
    {
        mesh_writer.AddCellData("Cell types", cell_types);
    }
    if (TissueConfig::Instance()->GetOutputCellAncestors())
    {
        mesh_writer.AddCellData("Ancestors", cell_ancestors);
    }
    if (TissueConfig::Instance()->GetOutputCellMutationStates())
    {
        mesh_writer.AddCellData("Mutation states", cell_mutation_states);
    }
    if (TissueConfig::Instance()->GetOutputCellAges())
    {
        mesh_writer.AddCellData("Ages", cell_ages);
    }
    if (TissueConfig::Instance()->GetOutputCellCyclePhases())
    {
        mesh_writer.AddCellData("Cycle phases", cell_cycle_phases);
    }
    if (TissueConfig::Instance()->GetOutputCellVolumes())
    {
        mesh_writer.AddCellData("Cell volumes", cell_volumes);
    }

    mesh_writer.WriteVtkUsingMesh(*(this->mpVoronoiTessellation), time.str());
    *(this->mpVtkMetaFile) << "        <DataSet timestep=\"";
    *(this->mpVtkMetaFile) << SimulationTime::Instance()->GetTimeStepsElapsed();
    *(this->mpVtkMetaFile) << "\" group=\"\" part=\"0\" file=\"results_";
    *(this->mpVtkMetaFile) << SimulationTime::Instance()->GetTimeStepsElapsed();
    *(this->mpVtkMetaFile) << ".vtu\"/>\n";
#endif //CHASTE_VTK
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////


template class MeshBasedTissueWithGhostNodes<1>;
template class MeshBasedTissueWithGhostNodes<2>;
template class MeshBasedTissueWithGhostNodes<3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MeshBasedTissueWithGhostNodes)
