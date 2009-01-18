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
#include "VertexBasedTissue.hpp"
#include "VertexMeshWriter.hpp"

template<unsigned DIM>
VertexBasedTissue<DIM>::VertexBasedTissue(VertexMesh<DIM, DIM>& rMesh,
                                          const std::vector<TissueCell>& rCells,
                                          bool deleteMesh,
                                          bool validate)
    : AbstractTissue<DIM>(rCells),
      mrMesh(rMesh),
      mDeleteMesh(deleteMesh)
{    
    // This must always be true
    assert( this->mCells.size() == mrMesh.GetNumElements() );

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
VertexElement<DIM, DIM>* VertexBasedTissue<DIM>::GetElementCorrespondingToCell(const TissueCell& rCell)
{
    return mrMesh.GetElement(rCell.GetLocationIndex());
}


template<unsigned DIM>
unsigned VertexBasedTissue<DIM>::GetNumElements()
{
    return mrMesh.GetNumElements();
}


template<unsigned DIM>
TissueCell* VertexBasedTissue<DIM>::AddCell(TissueCell& rNewCell, c_vector<double,DIM> newLocation)
{
//    // Get the element associated with this cell
//    unsigned element_index = GetElementCorrespondingToCell(rNewCell);    
//    VertexElement<DIM, DIM>* p_element = mrMesh.GetElement(element_index);
//
//    // Get the node indices owned by this element
//    std::set<unsigned> node_indices;
//    for (unsigned local_index=0; local_index<p_element->GetNumNodes(); local_index++)
//    {
//        node_indices.insert(p_element->GetNodeGlobalIndex(local_index));        
//    }
//
//    // Using this element's centroid and short axis, compute 
//    // the locations of two new nodes
//    std::vector<c_vector<double, DIM> > new_node_locations;
//    /// \todo Add method for computing new node locations here    
//    c_vector<double, DIM> newLocation1 = new_node_locations[0];
//    c_vector<double, DIM> newLocation2 = new_node_locations[1];
//    
//    // Create the two new nodes
//    Node<DIM>* p_new_node1 = new Node<DIM>(this->GetNumNodes(), newLocation1, false);
//    unsigned new_node1_index = AddNode(p_new_node1);
//    
//    Node<DIM>* p_new_node2 = new Node<DIM>(this->GetNumNodes(), newLocation2, false);
//    unsigned new_node1_index = AddNode(p_new_node1);
//    
//    /// \todo might the new nodes be on the boundary?
//
//    // Update the nodes owned by the existing element
//
//    // Add a new element to the mesh, which shared these two nodes
//    std::vector<Node<SPACE_DIM>*> nodes_in_new_element;
//    
//    /// \todo Write code to put nodes belonging to new element into the above vector
//    
//    VertexElement<DIM, DIM>* p_element = new VertexElement<DIM, DIM>(GetNumElements(), nodes_in_new_element);
//
//    // Associate the new cell with the element
//    newCell.SetLocationIndex(new_element_index);
//    this->mCells.push_back(rNewCell);
//
//    // Update location cell map
//    TissueCell *p_created_cell = &(this->mCells.back());
//    this->mLocationCellMap[new_element_index] = p_created_cell;
//
//    return p_created_cell;
//    
    return NULL; /// \todo put code for adding a cell here (see #852)
}


template<unsigned DIM>
unsigned VertexBasedTissue<DIM>::RemoveDeadCells()
{
    unsigned num_removed = 0;

    /// \todo put code for removing dead cells here (see #853)
    
    return num_removed;
}


template<unsigned DIM>
bool VertexBasedTissue<DIM>::IsCellAssociatedWithADeletedNode(TissueCell cell)
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
        this->mLocationCellMap.clear();
        
        for (std::list<TissueCell>::iterator cell_iter = this->mCells.begin();
             cell_iter != this->mCells.end();
             ++cell_iter)
        {
            unsigned old_elem_index = GetElementCorrespondingToCell(*cell_iter)->GetIndex();
            assert(!element_map.IsDeleted(old_elem_index));
            
            unsigned new_elem_index = element_map.GetNewIndex(old_elem_index);
            cell_iter->SetLocationIndex(new_elem_index);
            this->mLocationCellMap[new_elem_index] = &(*cell_iter);
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
        unsigned elem_index = GetElementCorrespondingToCell(*cell_iter)->GetIndex();
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
