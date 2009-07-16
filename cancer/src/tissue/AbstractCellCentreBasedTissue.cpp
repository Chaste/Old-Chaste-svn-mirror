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

#include "AbstractCellCentreBasedTissue.hpp"


template<unsigned DIM>
AbstractCellCentreBasedTissue<DIM>::AbstractCellCentreBasedTissue(const std::vector<TissueCell>& rCells,
                                                                  const std::vector<unsigned> locationIndices)
    : AbstractTissue<DIM>(rCells, locationIndices)
{
}


template<unsigned DIM>
AbstractCellCentreBasedTissue<DIM>::AbstractCellCentreBasedTissue()
    : AbstractTissue<DIM>()
{
}


template<unsigned DIM>
c_vector<double, DIM> AbstractCellCentreBasedTissue<DIM>::GetLocationOfCellCentre(TissueCell* pCell)
{
    return GetNodeCorrespondingToCell(pCell)->rGetLocation();
}


template<unsigned DIM>
Node<DIM>* AbstractCellCentreBasedTissue<DIM>::GetNodeCorrespondingToCell(TissueCell* pCell)
{
    return this->GetNode(this->mCellLocationMap[pCell]);
}


template<unsigned DIM>
TissueCell* AbstractCellCentreBasedTissue<DIM>::AddCell(TissueCell& rNewCell, c_vector<double,DIM> newLocation, TissueCell* pParentCell)
{
    // Create a new node
    Node<DIM> *p_new_node = new Node<DIM>(this->GetNumNodes(), newLocation, false);   // never on boundary
    unsigned new_node_index = AddNode(p_new_node); // use copy constructor so it doesn't matter that new_node goes out of scope

    // Update cells vector
    this->mCells.push_back(rNewCell);

    // Update mappings between cells and location indices
    TissueCell *p_created_cell = &(this->mCells.back());
    this->mLocationCellMap[new_node_index] = p_created_cell;
    this->mCellLocationMap[p_created_cell] = new_node_index;

    return p_created_cell;
}


template<unsigned DIM>
bool AbstractCellCentreBasedTissue<DIM>::IsCellAssociatedWithADeletedNode(TissueCell& rCell)
{
    return this->GetNode(this->mCellLocationMap[&rCell])->IsDeleted();
}


template<unsigned DIM>
void AbstractCellCentreBasedTissue<DIM>::UpdateNodeLocations(const std::vector< c_vector<double, DIM> >& rNodeForces, double dt)
{
    // Iterate over all nodes associated with real cells to update their positions
    for (typename AbstractTissue<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        // Get index of node associated with cell
        unsigned node_index = this->mCellLocationMap[&(*cell_iter)];

        // Get damping constant for node
        double damping_const = this->GetDampingConstant(node_index);

        // Get new node location
        c_vector<double, DIM> new_node_location = this->GetNode(node_index)->rGetLocation() + dt*rNodeForces[node_index]/damping_const;

        // Create ChastePoint for new node location
        ChastePoint<DIM> new_point(new_node_location);

        // Move the node
        this->SetNode(node_index, new_point);
    }
}


template<unsigned DIM>
double AbstractCellCentreBasedTissue<DIM>::GetDampingConstant(unsigned nodeIndex)
{
    double damping_multiplier = 1.0;

    if (   (this->rGetCellUsingLocationIndex(nodeIndex).GetMutationState() != HEALTHY)
        && (this->rGetCellUsingLocationIndex(nodeIndex).GetMutationState() != APC_ONE_HIT) )
    {
        return TissueConfig::Instance()->GetDampingConstantMutant()*damping_multiplier;
    }
    else
    {
        return TissueConfig::Instance()->GetDampingConstantNormal()*damping_multiplier;
    }
}


template<unsigned DIM>
void AbstractCellCentreBasedTissue<DIM>::WriteResultsToFiles()
{
    AbstractTissue<DIM>::WriteResultsToFiles();

    std::vector<unsigned> cell_type_counter, cell_mutation_state_counter, cell_cycle_phase_counter;

    this->WriteTimeAndNodeResultsToFiles(cell_type_counter,
                                         cell_mutation_state_counter,
                                         cell_cycle_phase_counter);

	/*
	 * Note that we cannot use a mesh iterator below, as the child class
	 * NodeBasedTissue does not contain a mesh.
	 */
    for (unsigned node_index=0; node_index<this->GetNumNodes(); node_index++)
    {
        if ( !(this->GetNode(node_index)->IsDeleted()) )
        {
            this->GenerateCellResults(node_index,
                                      cell_type_counter,
                                      cell_mutation_state_counter,
                                      cell_cycle_phase_counter);
        }
    }

    this->WriteCellResultsToFiles(cell_type_counter,
                                  cell_mutation_state_counter,
                                  cell_cycle_phase_counter);
}


/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class AbstractCellCentreBasedTissue<1>;
template class AbstractCellCentreBasedTissue<2>;
template class AbstractCellCentreBasedTissue<3>;
