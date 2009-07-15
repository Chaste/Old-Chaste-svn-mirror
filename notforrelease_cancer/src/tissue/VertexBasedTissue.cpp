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

#include "VertexBasedTissue.hpp"
#include "VertexMeshWriter.hpp"

template<unsigned DIM>
VertexBasedTissue<DIM>::VertexBasedTissue(VertexMesh<DIM, DIM>& rMesh,
                                          const std::vector<TissueCell>& rCells,
                                          bool deleteMesh,
                                          bool validate,
                                          const std::vector<unsigned> locationIndices)
    : AbstractTissue<DIM>(rCells, locationIndices),
      mrMesh(rMesh),
      mDeleteMesh(deleteMesh)
{
    this->mTissueContainsMesh = true;

    if (validate)
    {
        Validate();
    }
}


template<unsigned DIM>
VertexBasedTissue<DIM>::VertexBasedTissue(VertexMesh<DIM, DIM>& rMesh)
             : mrMesh(rMesh)
{
    this->mTissueContainsMesh = true;
    mDeleteMesh = true;
}


template<unsigned DIM>
VertexBasedTissue<DIM>::~VertexBasedTissue()
{
    if (mDeleteMesh)
    {
        delete &mrMesh;
    }
}


template<unsigned DIM>
double VertexBasedTissue<DIM>::GetDampingConstant(unsigned nodeIndex)
{
    // Take the average of the cells containing this vertex
    double average_damping_constant = 0.0;

    std::set<unsigned> containing_elements = GetNode(nodeIndex)->rGetContainingElementIndices();
    unsigned num_containing_elements = containing_elements.size();

    for (std::set<unsigned>::iterator iter=containing_elements.begin();
         iter!=containing_elements.end();
         ++iter)
    {
        if (   (this->rGetCellUsingLocationIndex(*iter).GetMutationState() != HEALTHY)
            && (this->rGetCellUsingLocationIndex(*iter).GetMutationState() != APC_ONE_HIT) )
        {
            average_damping_constant += TissueConfig::Instance()->GetDampingConstantMutant()/((double) num_containing_elements);
        }
        else
        {
            average_damping_constant += TissueConfig::Instance()->GetDampingConstantNormal()/((double) num_containing_elements);
        }
    }

    return average_damping_constant;
}


template<unsigned DIM>
double VertexBasedTissue<DIM>::GetAdhesionParameter(Node<DIM>* pNodeA, Node<DIM>* pNodeB)
{
    double adhesion_parameter;

    // Find the indices of the elements owned by each node
    std::set<unsigned> elements_containing_nodeA = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> elements_containing_nodeB = pNodeB->rGetContainingElementIndices();

    // Find common elements
    std::set<unsigned> shared_elements;
    std::set_intersection(elements_containing_nodeA.begin(),
                          elements_containing_nodeA.end(),
                          elements_containing_nodeB.begin(),
                          elements_containing_nodeB.end(),
                          std::inserter(shared_elements, shared_elements.begin()));

    // Check that the nodes have a common edge
    assert(shared_elements.size() > 0);

    // If the edge corresponds to a single element, then the cell is on the boundary
    if (shared_elements.size() == 1)
    {
        adhesion_parameter = TissueConfig::Instance()->GetCellBoundaryAdhesionEnergyParameter();
    }
    else
    {
        adhesion_parameter = TissueConfig::Instance()->GetCellCellAdhesionEnergyParameter();
    }
    return adhesion_parameter;
}


template<unsigned DIM>
VertexMesh<DIM, DIM>& VertexBasedTissue<DIM>::rGetMesh()
{
    return mrMesh;
}


template<unsigned DIM>
const VertexMesh<DIM, DIM>& VertexBasedTissue<DIM>::rGetMesh() const
{
    return mrMesh;
}


template<unsigned DIM>
VertexElement<DIM, DIM>* VertexBasedTissue<DIM>::GetElement(unsigned elementIndex)
{
    return mrMesh.GetElement(elementIndex);
}


template<unsigned DIM>
unsigned VertexBasedTissue<DIM>::GetNumNodes()
{
    return mrMesh.GetNumNodes();
}


template<unsigned DIM>
c_vector<double, DIM> VertexBasedTissue<DIM>::GetLocationOfCellCentre(TissueCell* pCell)
{
    // Get location index corresponding to this cell
    unsigned location_index = this->GetLocationIndexUsingCell(pCell);
    return mrMesh.GetCentroidOfElement(location_index);
}


template<unsigned DIM>
Node<DIM>* VertexBasedTissue<DIM>::GetNode(unsigned index)
{
    return mrMesh.GetNode(index);
}


template<unsigned DIM>
unsigned VertexBasedTissue<DIM>::AddNode(Node<DIM>* pNewNode)
{
    return mrMesh.AddNode(pNewNode);
}


template<unsigned DIM>
void VertexBasedTissue<DIM>::SetNode(unsigned nodeIndex, ChastePoint<DIM>& rNewLocation)
{
    mrMesh.SetNode(nodeIndex, rNewLocation);
}


template<unsigned DIM>
VertexElement<DIM, DIM>* VertexBasedTissue<DIM>::GetElementCorrespondingToCell(TissueCell* pCell)
{
    return mrMesh.GetElement(this->mCellLocationMap[pCell]);
}


template<unsigned DIM>
unsigned VertexBasedTissue<DIM>::GetNumElements()
{
    return mrMesh.GetNumElements();
}


template<unsigned DIM>
TissueCell* VertexBasedTissue<DIM>::AddCell(TissueCell& rNewCell, c_vector<double,DIM> newLocation, TissueCell* pParentCell)
{
    /// \todo newLocation is redundant for vertex-based tissues (#852)

    // Get the element associated with this cell
    VertexElement<DIM, DIM> *p_element = GetElementCorrespondingToCell(pParentCell);

    // Divide the element
    unsigned new_element_index; 
    if (pParentCell->GetCellType() == STEM)
    {
        // Divide this element horizontally
        c_vector<double, 2> axis_of_division;
        axis_of_division(0)=1.0;
        axis_of_division(1)=0.0;
        
        new_element_index = mrMesh.DivideElement(p_element,axis_of_division);
    }
    else
    {
        // Divide this element along the short axis
        new_element_index = mrMesh.DivideElement(p_element);
    }

    // Associate the new cell with the element
    this->mCells.push_back(rNewCell);

    // Update location cell map
    TissueCell *p_created_cell = &(this->mCells.back());
    this->mLocationCellMap[new_element_index] = p_created_cell;
    this->mCellLocationMap[p_created_cell] = new_element_index;

    return p_created_cell;
}


template<unsigned DIM>
unsigned VertexBasedTissue<DIM>::RemoveDeadCells()
{
    //std::cout << "time, apop = " << SimulationTime::Instance()->GetTime() << " " << this->rGetCellUsingLocationIndex(18).HasApoptosisBegun() << "\n";
        
    unsigned num_removed = 0;

    for (std::list<TissueCell>::iterator it = this->mCells.begin();
         it != this->mCells.end();
         ++it)
    {
        if (it->IsDead())
        {
            // Remove the element from the mesh
            num_removed++;
            mrMesh.DeleteElementPriorToReMesh(this->mCellLocationMap[&(*it)]);
            it = this->mCells.erase(it);
            --it;
        }
    }
    return num_removed;
}


template<unsigned DIM>
void VertexBasedTissue<DIM>::UpdateNodeLocations(const std::vector< c_vector<double, DIM> >& rNodeForces, double dt)
{
    // Iterate over all nodes associated with real cells to update their positions
    for (unsigned node_index=0; node_index<GetNumNodes(); node_index++)
    {
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
bool VertexBasedTissue<DIM>::IsCellAssociatedWithADeletedNode(TissueCell& rCell)
{
    return false;
}


template<unsigned DIM>
void VertexBasedTissue<DIM>::Update(bool hasHadBirthsOrDeaths)
{
    VertexElementMap element_map(mrMesh.GetNumAllElements());

    mrMesh.ReMesh(element_map);

    if (!element_map.IsIdentityMap())
    {
        // Fix up the mappings between TissueCells and VertexElements
        std::map<TissueCell*, unsigned> old_map = this->mCellLocationMap;

        this->mCellLocationMap.clear();
        this->mLocationCellMap.clear();

        for (std::list<TissueCell>::iterator cell_iter = this->mCells.begin();
             cell_iter != this->mCells.end();
             ++cell_iter)
        {
            // This shouldn't ever happen, as the cell vector only contains living cells
            unsigned old_elem_index = old_map[&(*cell_iter)];
            assert(!element_map.IsDeleted(old_elem_index));
            unsigned new_elem_index = element_map.GetNewIndex(old_elem_index);

            this->mLocationCellMap[new_elem_index] = &(*cell_iter);
            this->mCellLocationMap[&(*cell_iter)] = new_elem_index;
        }
    }
    // Check that each VertexElement has a TissueCell associated with it in the updated tissue
    Validate();
}


template<unsigned DIM>
void VertexBasedTissue<DIM>::Validate()
{
    std::vector<bool> validated_element = std::vector<bool>(this->GetNumElements(), false);

    for (typename AbstractTissue<DIM>::Iterator cell_iter=this->Begin();
         cell_iter!=this->End();
         ++cell_iter)
    {
        unsigned elem_index = GetLocationIndexUsingCell(&(*cell_iter));
        validated_element[elem_index] = true;
    }

    for (unsigned i=0; i<validated_element.size(); i++)
    {
        if (!validated_element[i])
        {
            std::stringstream ss;
            ss << "Element " << i << " does not appear to have a cell associated with it";
            EXCEPTION(ss.str());
        }
    }
}


template<unsigned DIM>
double VertexBasedTissue<DIM>::GetTargetAreaOfCell(const TissueCell& rCell)
{
    // Get target area A of a healthy cell in S, G2 or M phase
    double cell_target_area = TissueConfig::Instance()->GetMatureCellTargetArea();

    double cell_age = rCell.GetAge();
    double g1_duration = rCell.GetCellCycleModel()->GetG1Duration();
        
    // If differentiated then g1_duration is infinite
    if (g1_duration == DBL_MAX) // dont use magic number, compare to DBL_MAX
    {
        // This is just for fixed cell cycle models, need to work out how to find the g1 duration
        g1_duration = TissueConfig::Instance()->GetTransitCellG1Duration();
    }    

    if (rCell.GetCellType()==APOPTOTIC)
    {
        // Age of cell when apoptosis begins
        //double cell_age_at_death = cell_age - TissueConfig::Instance()->GetApoptosisTime() + rCell.TimeUntilDeath(); 
        
        if (rCell.GetStartOfApoptosisTime() - rCell.GetBirthTime() < g1_duration)
        {
            cell_target_area *= 0.5*(1 + (rCell.GetStartOfApoptosisTime() - rCell.GetBirthTime())/g1_duration);
        }

        // The target area of an apoptotic cell decreases linearly to zero (and past it negative)        
        cell_target_area = cell_target_area - cell_target_area/(TissueConfig::Instance()->GetApoptosisTime())*(SimulationTime::Instance()->GetTime()-rCell.GetStartOfApoptosisTime());
    
        // Don't allow a negative target area
        if (cell_target_area < 0)
        {
            cell_target_area = 0;
        }
    }
    else 
    {
        // In the case of a proliferating cell, the target area increases 
        // linearly from A/2 to A over the course of the G1 phase 
        if (cell_age < g1_duration)
        {
            cell_target_area *= 0.5*(1 + cell_age/g1_duration);
        }
    }

    return cell_target_area;
}

template<unsigned DIM>
void VertexBasedTissue<DIM>::WriteResultsToFiles(bool outputCellMutationStates,
                                                 bool outputCellTypes,
                                                 bool outputCellVariables,
                                                 bool outputCellCyclePhases,
                                                 bool outputCellAncestors,
                                                 bool outputCellAges)
{
    std::vector<unsigned> cell_type_counter, cell_mutation_state_counter, cell_cycle_phase_counter;

    this->WriteTimeAndNodeResultsToFiles(outputCellMutationStates,
                                         outputCellTypes,
                                         outputCellVariables,
                                         outputCellCyclePhases,
                                         outputCellAncestors,
                                         outputCellAges,
                                         cell_type_counter,
                                         cell_mutation_state_counter,
                                         cell_cycle_phase_counter);

    // Write element data to file
    *mpElementFile << SimulationTime::Instance()->GetTime() << "\t";
    for (typename VertexMesh<DIM,DIM>::VertexElementIterator iter = mrMesh.GetElementIteratorBegin();
             iter != mrMesh.GetElementIteratorEnd();
             ++iter)
    {
        if ( !(iter->IsDeleted()) )
        {
            // First write the number of Nodes belonging to this VertexElement
            *mpElementFile << iter->GetNumNodes() << " ";

            // Then write the global index of each Node in this element
            unsigned num_nodes_in_this_element = iter->GetNumNodes();
            for (unsigned i=0; i<num_nodes_in_this_element; i++)
            {
                *mpElementFile << iter->GetNodeGlobalIndex(i) << " ";
            }
        }
    }
    *mpElementFile << "\n";

    for (typename VertexMesh<DIM,DIM>::VertexElementIterator iter = mrMesh.GetElementIteratorBegin();
             iter != mrMesh.GetElementIteratorEnd();
             ++iter)
    {
        if (!(iter->IsDeleted()))
        {
            this->GenerateCellResults(iter->GetIndex(),
                                      outputCellMutationStates,
                                      outputCellTypes,
                                      outputCellVariables,
                                      outputCellCyclePhases,
                                      outputCellAncestors,
                                      outputCellAges,
                                      cell_type_counter,
                                      cell_mutation_state_counter,
                                      cell_cycle_phase_counter);
        }
    }

    this->WriteCellResultsToFiles(outputCellMutationStates,
                                  outputCellTypes,
                                  outputCellVariables,
                                  outputCellCyclePhases,
                                  outputCellAncestors,
                                  outputCellAges,
                                  cell_type_counter,
                                  cell_mutation_state_counter,
                                  cell_cycle_phase_counter);
}


template<unsigned DIM>
void VertexBasedTissue<DIM>::WriteMeshToFile(const std::string& rArchiveDirectory, const std::string& rMeshFileName)
{
    // The false is so the directory isn't cleaned
    VertexMeshWriter<DIM, DIM> mesh_writer(rArchiveDirectory, rMeshFileName, false);

    mesh_writer.WriteFilesUsingMesh(mrMesh);
}


template<unsigned DIM>
void VertexBasedTissue<DIM>::CreateOutputFiles(const std::string& rDirectory,
                                               bool rCleanOutputDirectory,
                                               bool outputCellMutationStates,
                                               bool outputCellTypes,
                                               bool outputCellVariables,
                                               bool outputCellCyclePhases,
                                               bool outputCellAncestors,
                                               bool outputCellAges)
{
    AbstractTissue<DIM>::CreateOutputFiles(rDirectory,
                                           rCleanOutputDirectory,
                                           outputCellMutationStates,
                                           outputCellTypes,
                                           outputCellVariables,
                                           outputCellCyclePhases,
                                           outputCellAncestors,
                                           outputCellAges);

    OutputFileHandler output_file_handler(rDirectory, rCleanOutputDirectory);
    mpElementFile = output_file_handler.OpenOutputFile("results.vizelements");
}


template<unsigned DIM>
void VertexBasedTissue<DIM>::CloseOutputFiles(bool outputCellMutationStates,
                                              bool outputCellTypes,
                                              bool outputCellVariables,
                                              bool outputCellCyclePhases,
                                              bool outputCellAncestors,
                                              bool outputCellAges)
{
    AbstractTissue<DIM>::CloseOutputFiles(outputCellMutationStates,
                                          outputCellTypes,
                                          outputCellVariables,
                                          outputCellCyclePhases,
                                          outputCellAncestors,
                                          outputCellAges);
    mpElementFile->close();
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////


template class VertexBasedTissue<1>;
template class VertexBasedTissue<2>;
template class VertexBasedTissue<3>;
