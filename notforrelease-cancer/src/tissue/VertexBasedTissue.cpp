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
void VertexBasedTissue<DIM>::SetNode(unsigned index, ChastePoint<DIM>& rNewLocation)
{
    return mrMesh.SetNode(index, rNewLocation);
}


template<unsigned DIM>
unsigned VertexBasedTissue<DIM>::GetNumElements()
{
    return mrMesh.GetNumElements();
}


template<unsigned DIM>
void VertexBasedTissue<DIM>::MoveCell(typename AbstractTissue<DIM>::Iterator iter, ChastePoint<DIM>& rNewLocation)
{
    /// \todo put code for moving a cell here
}


template<unsigned DIM>
TissueCell* VertexBasedTissue<DIM>::AddCell(TissueCell cell, c_vector<double,DIM> newLocation)
{
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
        
        for (std::list<TissueCell>::iterator iter = this->mCells.begin();
             iter != this->mCells.end();
             ++iter)
        {
            unsigned old_elem_index = iter->GetLocationIndex();
            assert(!element_map.IsDeleted(old_elem_index));
            
            unsigned new_elem_index = element_map.GetNewIndex(old_elem_index);
            iter->SetLocationIndex(new_elem_index);
            this->mLocationCellMap[new_elem_index] = &(*iter);
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
        unsigned elem_index = cell_iter->GetLocationIndex();
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


template<unsigned DIM>
std::string VertexBasedTissue<DIM>::meshPathname = "";


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////


template class VertexBasedTissue<1>;
template class VertexBasedTissue<2>;
template class VertexBasedTissue<3>;
