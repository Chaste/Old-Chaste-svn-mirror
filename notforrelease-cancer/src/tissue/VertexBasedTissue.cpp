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
    // This must always be true
    if (this->mCells.size() != mrMesh.GetNumElements())
    {
        std::stringstream ss;
        ss << "The number of cells does not match the number of elements in the mesh";
        EXCEPTION(ss.str());
    }

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
    /// \todo Implement variable damping constants for a vertex-based tissue (see #865)
    return CancerParameters::Instance()->GetDampingConstantNormal();
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
Node<DIM>* VertexBasedTissue<DIM>::GetNode(unsigned index)
{
    return mrMesh.GetNode(index);
}


template<unsigned DIM>
unsigned VertexBasedTissue<DIM>::AddNode(Node<DIM> *pNewNode)
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
    VertexElement<DIM, DIM>* p_element = GetElementCorrespondingToCell(pParentCell);

    // Divide this element along the short axis
    unsigned new_element_index = mrMesh.DivideElement(p_element);

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
    unsigned num_removed = 0;

    /// \todo put code for removing dead cells here (see #853)

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
void VertexBasedTissue<DIM>::Update()
{
    /// \todo Thought about creating an ElementMap class, but it looks like
    //        we can just hijack NodeMap for our purposes... in fact, what *is*
    //        specific to Nodes in NodeMap?? (see #827)
    NodeMap element_map(mrMesh.GetNumElements());
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
        unsigned elem_index = GetElementCorrespondingToCell(&(*cell_iter))->GetIndex();
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
void VertexBasedTissue<DIM>::WriteResultsToFiles(bool outputCellMutationStates,
                                                 bool outputCellTypes,
                                                 bool outputCellVariables,
                                                 bool outputCellCyclePhases,
                                                 bool outputCellAncestors)
{
    std::vector<unsigned> cell_type_counter, cell_mutation_state_counter, cell_cycle_phase_counter;

    this->WriteTimeAndNodeResultsToFiles(outputCellMutationStates,
                                         outputCellTypes,
                                         outputCellVariables,
                                         outputCellCyclePhases,
                                         outputCellAncestors,
                                         cell_type_counter,
                                         cell_mutation_state_counter,
                                         cell_cycle_phase_counter);

    // Write element data to file
    *mpElementFile << SimulationTime::Instance()->GetTime() << "\t";
    for (unsigned elem_index=0; elem_index<GetNumElements(); elem_index++)
    {
        if ( !(mrMesh.GetElement(elem_index)->IsDeleted()) )
        {
            // First write the number of Nodes belonging to this VertexElement
            *mpElementFile << mrMesh.GetElement(elem_index)->GetNumNodes() << " ";

            // Then write the global index of each Node in this element
            unsigned num_nodes_in_this_element = mrMesh.GetElement(elem_index)->GetNumNodes();
            for (unsigned i=0; i<num_nodes_in_this_element; i++)
            {
                *mpElementFile << mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(i) << " ";
            }
        }
    }
    *mpElementFile << "\n";

    for (unsigned elem_index=0; elem_index<GetNumElements(); elem_index++)
    {
        if (!(mrMesh.GetElement(elem_index)->IsDeleted()))
        {
            this->GenerateCellResults(elem_index,
                                      outputCellMutationStates,
                                      outputCellTypes,
                                      outputCellVariables,
                                      outputCellCyclePhases,
                                      outputCellAncestors,
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
                                  cell_type_counter,
                                  cell_mutation_state_counter,
                                  cell_cycle_phase_counter);
}


template<unsigned DIM>
void VertexBasedTissue<DIM>::WriteMeshToFile(const std::string &rArchiveDirectory, const std::string &rMeshFileName)
{
    // The false is so the directory isn't cleaned
    VertexMeshWriter<DIM, DIM> mesh_writer(rArchiveDirectory, rMeshFileName, false);

    mesh_writer.WriteFilesUsingMesh(mrMesh);
}


template<unsigned DIM>
void VertexBasedTissue<DIM>::CreateOutputFiles(const std::string &rDirectory,
                                               bool rCleanOutputDirectory,
                                               bool outputCellMutationStates,
                                               bool outputCellTypes,
                                               bool outputCellVariables,
                                               bool outputCellCyclePhases,
                                               bool outputCellAncestors)
{
    AbstractTissue<DIM>::CreateOutputFiles(rDirectory,
                                           rCleanOutputDirectory,
                                           outputCellMutationStates,
                                           outputCellTypes,
                                           outputCellVariables,
                                           outputCellCyclePhases,
                                           outputCellAncestors);

    OutputFileHandler output_file_handler(rDirectory, rCleanOutputDirectory);
    mpElementFile = output_file_handler.OpenOutputFile("results.vizelements");
}


template<unsigned DIM>
void VertexBasedTissue<DIM>::CloseOutputFiles(bool outputCellMutationStates,
                                              bool outputCellTypes,
                                              bool outputCellVariables,
                                              bool outputCellCyclePhases,
                                              bool outputCellAncestors)
{
    AbstractTissue<DIM>::CloseOutputFiles(outputCellMutationStates,
                                          outputCellTypes,
                                          outputCellVariables,
                                          outputCellCyclePhases,
                                          outputCellAncestors);
    mpElementFile->close();
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////


template class VertexBasedTissue<1>;
template class VertexBasedTissue<2>;
template class VertexBasedTissue<3>;
