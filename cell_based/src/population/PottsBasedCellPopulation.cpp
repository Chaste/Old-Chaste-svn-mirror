/*

Copyright (c) 2005-2012, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "PottsBasedCellPopulation.hpp"
#include "RandomNumberGenerator.hpp"
#include "Warnings.hpp"

// Needed to convert mesh in order to write nodes to VTK (visualize as glyphs)
#include "VtkMeshWriter.hpp"
#include "NodesOnlyMesh.hpp"
#include "Exception.hpp"

template<unsigned DIM>
void PottsBasedCellPopulation<DIM>::Validate()
{
    // Check each element has only one cell associated with it
    std::vector<unsigned> validated_element = std::vector<unsigned>(this->GetNumElements(), 0);

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        unsigned elem_index = this->GetLocationIndexUsingCell(*cell_iter);
        validated_element[elem_index]++;
    }

    for (unsigned i=0; i<validated_element.size(); i++)
    {
        if (validated_element[i] == 0)
        {
            EXCEPTION("Element " << i << " does not appear to have a cell associated with it");
        }

        if (validated_element[i] > 1)
        {
            EXCEPTION("Element " << i << " appears to have " << validated_element[i] << " cells associated with it");
        }
    }
}

template<unsigned DIM>
PottsBasedCellPopulation<DIM>::PottsBasedCellPopulation(PottsMesh<DIM>& rMesh,
                                                        std::vector<CellPtr>& rCells,
                                                        bool deleteMesh,
                                                        bool validate,
                                                        const std::vector<unsigned> locationIndices)
    : AbstractOnLatticeCellPopulation<DIM>(rMesh, rCells, locationIndices, deleteMesh),
      mpElementTessellation(NULL),
      mpMutableMesh(NULL),
      mTemperature(0.1),
      mNumSweepsPerTimestep(1)
{
	mpPottsMesh = static_cast<PottsMesh<DIM>* >(&(this->mrMesh));
    // Check each element has only one cell associated with it
    if (validate)
    {
        Validate();
    }
}

template<unsigned DIM>
PottsBasedCellPopulation<DIM>::PottsBasedCellPopulation(PottsMesh<DIM>& rMesh)
    : AbstractOnLatticeCellPopulation<DIM>(rMesh),
      mpElementTessellation(NULL),
      mpMutableMesh(NULL),
      mTemperature(0.1),
      mNumSweepsPerTimestep(1)
{
	mpPottsMesh = static_cast<PottsMesh<DIM>* >(&(this->mrMesh));
}

template<unsigned DIM>
PottsBasedCellPopulation<DIM>::~PottsBasedCellPopulation()
{
    delete mpElementTessellation;

    delete mpMutableMesh;

    if (this->mDeleteMesh)
    {
        delete &this->mrMesh;
    }
}

template<unsigned DIM>
PottsMesh<DIM>& PottsBasedCellPopulation<DIM>::rGetMesh()
{
    return *mpPottsMesh;
}

template<unsigned DIM>
const PottsMesh<DIM>& PottsBasedCellPopulation<DIM>::rGetMesh() const
{
    return *mpPottsMesh;
}

template<unsigned DIM>
PottsElement<DIM>* PottsBasedCellPopulation<DIM>::GetElement(unsigned elementIndex)
{
    return mpPottsMesh->GetElement(elementIndex);
}

template<unsigned DIM>
unsigned PottsBasedCellPopulation<DIM>::GetNumElements()
{
    return mpPottsMesh->GetNumElements();
}

template<unsigned DIM>
Node<DIM>* PottsBasedCellPopulation<DIM>::GetNode(unsigned index)
{
    return this->mrMesh.GetNode(index);
}

template<unsigned DIM>
unsigned PottsBasedCellPopulation<DIM>::GetNumNodes()
{
    return this->mrMesh.GetNumNodes();
}

template<unsigned DIM>
c_vector<double, DIM> PottsBasedCellPopulation<DIM>::GetLocationOfCellCentre(CellPtr pCell)
{
    return mpPottsMesh->GetCentroidOfElement(this->GetLocationIndexUsingCell(pCell));
}

template<unsigned DIM>
PottsElement<DIM>* PottsBasedCellPopulation<DIM>::GetElementCorrespondingToCell(CellPtr pCell)
{
    return mpPottsMesh->GetElement(this->GetLocationIndexUsingCell(pCell));
}

template<unsigned DIM>
CellPtr PottsBasedCellPopulation<DIM>::AddCell(CellPtr pNewCell, const c_vector<double,DIM>& rCellDivisionVector, CellPtr pParentCell)
{
    // Get the element associated with this cell
    PottsElement<DIM>* p_element = GetElementCorrespondingToCell(pParentCell);

    // Divide the element
    unsigned new_element_index = mpPottsMesh->DivideElement(p_element, true); // new element will be below the existing element

    // Associate the new cell with the element
    this->mCells.push_back(pNewCell);

    // Update location cell map
    CellPtr p_created_cell = this->mCells.back();
    this->SetCellUsingLocationIndex(new_element_index,p_created_cell);
    return p_created_cell;
}

template<unsigned DIM>
unsigned PottsBasedCellPopulation<DIM>::RemoveDeadCells()
{
    unsigned num_removed = 0;

    for (std::list<CellPtr>::iterator it = this->mCells.begin();
         it != this->mCells.end();
         ++it)
    {
        if ((*it)->IsDead())
        {
            // Remove the element from the mesh
            num_removed++;
            mpPottsMesh->DeleteElement(this->GetLocationIndexUsingCell((*it)));
            it = this->mCells.erase(it);
            --it;
        }
    }
    return num_removed;
}
template<unsigned DIM>
void PottsBasedCellPopulation<DIM>::UpdateCellLocations(double dt)
{
    /*
     * This method implements a Monte Carlo method to update the cell population.
     * We sample randomly from all nodes in the mesh. Once we have selected a target
     * node we randomly select a neighbour. The Hamiltonian is evaluated in the
     * current configuration (H_0) and with the target node added to the same
     * element as the neighbour (H_1). Based on the vale of deltaH = H_1 - H_0,
     * the switch is either made or not.
     *
     * For each time step (i.e. each time this method is called) we sample
     * mrMesh.GetNumNodes() nodes. This is known as a Monte Carlo Step (MCS).
     */

    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    unsigned num_nodes = this->mrMesh.GetNumNodes();

    // Randomly permute mUpdateRuleCollection if specified
    if (this->mIterateRandomlyOverUpdateRuleCollection)
    {
        // Randomly permute mUpdateRuleCollection
    	p_gen->Shuffle(mUpdateRuleCollection);
    }

    for (unsigned i=0; i<num_nodes*mNumSweepsPerTimestep; i++)
    {
        unsigned node_index;

        if (this->mUpdateNodesInRandomOrder)
        {
            node_index = p_gen->randMod(num_nodes);
        }
        else
        {
            // Loop over nodes in index order.
            node_index = i%num_nodes;
        }

        Node<DIM>* p_node = this->mrMesh.GetNode(node_index);

        // Each node in the mesh must be in at most one element
        assert(p_node->GetNumContainingElements() <= 1);

        // Find a random available neighbouring node to overwrite current site
        std::set<unsigned> neighbouring_node_indices = mpPottsMesh->GetMooreNeighbouringNodeIndices(node_index);
        unsigned neighbour_location_index;

        if (!neighbouring_node_indices.empty())
        {
            unsigned num_neighbours = neighbouring_node_indices.size();
            unsigned chosen_neighbour = p_gen->randMod(num_neighbours);

            std::set<unsigned>::iterator neighbour_iter = neighbouring_node_indices.begin();
            for (unsigned i=0; i<chosen_neighbour; i++)
            {
                neighbour_iter++;
            }

            neighbour_location_index = *neighbour_iter;
        }
        else
        {
            // Each node in the mesh must have at least one neighbour
            NEVER_REACHED;
        }

        std::set<unsigned> containing_elements = p_node->rGetContainingElementIndices();
        std::set<unsigned> neighbour_containing_elements = GetNode(neighbour_location_index)->rGetContainingElementIndices();
        // Only calculate Hamiltonian and update elements if the nodes are from different elements, or one is from the medium
        if (  (  *containing_elements.begin() != *neighbour_containing_elements.begin() && !containing_elements.empty() && !neighbour_containing_elements.empty() )
                || ( !containing_elements.empty() && neighbour_containing_elements.empty() )
                || ( containing_elements.empty() && !neighbour_containing_elements.empty() ) )
        {
            double delta_H = 0.0; // This is H_1-H_0.

            // Now add contributions to the Hamiltonian from each AbstractPottsUpdateRule
            for (typename std::vector<boost::shared_ptr<AbstractPottsUpdateRule<DIM> > >::iterator iter = mUpdateRuleCollection.begin();
                 iter != mUpdateRuleCollection.end();
                 ++iter)
            {
                delta_H += (*iter)->EvaluateHamiltonianContribution(neighbour_location_index, p_node->GetIndex(), *this);
            }

            // Generate a uniform random number to do the random motion
            double random_number = p_gen->ranf();
            double p = exp(-delta_H/mTemperature);
            if (delta_H <= 0 || random_number < p)
            {
                // Do swap

                // Remove the current node from any elements containing it (there should be at most one such element)
                for (std::set<unsigned>::iterator iter = containing_elements.begin();
                     iter != containing_elements.end();
                     ++iter)
                {
                    GetElement(*iter)->DeleteNode(GetElement(*iter)->GetNodeLocalIndex(node_index));

                    ///\todo If this causes the element to have no nodes then flag the element and cell to be deleted
                }

                // Next add the current node to any elements containing the neighbouring node (there should be at most one such element)
                for (std::set<unsigned>::iterator iter = neighbour_containing_elements.begin();
                     iter != neighbour_containing_elements.end();
                     ++iter)
                {
                    GetElement(*iter)->AddNode(this->mrMesh.GetNode(node_index));
                }
            }
        }
    }
}

template<unsigned DIM>
bool PottsBasedCellPopulation<DIM>::IsCellAssociatedWithADeletedLocation(CellPtr pCell)
{
    return GetElementCorrespondingToCell(pCell)->IsDeleted();
}

template<unsigned DIM>
void PottsBasedCellPopulation<DIM>::Update(bool hasHadBirthsOrDeaths)
{
}

template<unsigned DIM>
void PottsBasedCellPopulation<DIM>::CreateOutputFiles(const std::string& rDirectory, bool cleanOutputDirectory)
{
    AbstractCellPopulation<DIM>::CreateOutputFiles(rDirectory, cleanOutputDirectory);

    OutputFileHandler output_file_handler(rDirectory, cleanOutputDirectory);
    mpVizElementsFile = output_file_handler.OpenOutputFile("results.vizelements");
}

template<unsigned DIM>
void PottsBasedCellPopulation<DIM>::CloseOutputFiles()
{
    AbstractCellPopulation<DIM>::CloseOutputFiles();
    mpVizElementsFile->close();
}

template<unsigned DIM>
void PottsBasedCellPopulation<DIM>::WriteResultsToFiles()
{
    AbstractCellPopulation<DIM>::WriteResultsToFiles();

    CreateElementTessellation(); // To be used to output to the visualizer

    SimulationTime* p_time = SimulationTime::Instance();

    // Write element data to file
    *mpVizElementsFile << p_time->GetTime() << "\t";

    // Loop over cells and find associated elements so in the same order as the cells in output files
    for (std::list<CellPtr>::iterator cell_iter = this->mCells.begin();
         cell_iter != this->mCells.end();
         ++cell_iter)
    {
        unsigned elem_index = this->GetLocationIndexUsingCell(*cell_iter);

        // Hack that covers the case where the element is associated with a cell that has just been killed (#1129)
        bool elem_corresponds_to_dead_cell = false;

        if (this->IsCellAttachedToLocationIndex(elem_index))
        {
            elem_corresponds_to_dead_cell = this->GetCellUsingLocationIndex(elem_index)->IsDead();
        }

        // Write node data to file
        if (!(GetElement(elem_index)->IsDeleted()) && !elem_corresponds_to_dead_cell)
        {
            PottsElement<DIM>* p_element = mpPottsMesh->GetElement(elem_index);

            unsigned num_nodes_in_element = p_element->GetNumNodes();

            // First write the number of Nodes belonging to this PottsElement
            *mpVizElementsFile << num_nodes_in_element << " ";

            // Then write the global index of each Node in this element
            for (unsigned i=0; i<num_nodes_in_element; i++)
            {
                *mpVizElementsFile << p_element->GetNodeGlobalIndex(i) << " ";
            }
        }
    }
    *mpVizElementsFile << "\n";
}

template<unsigned DIM>
double PottsBasedCellPopulation<DIM>::GetVolumeOfCell(CellPtr pCell)
{
    // Get element index corresponding to this cell
    unsigned elem_index = this->GetLocationIndexUsingCell(pCell);

    // Get volume of this element in the Potts mesh
    double cell_volume = mpPottsMesh->GetVolumeOfElement(elem_index);

    return cell_volume;
}

template<unsigned DIM>
void PottsBasedCellPopulation<DIM>::WriteCellVolumeResultsToFile()
{
    // Write time to file
    *(this->mpCellVolumesFile) << SimulationTime::Instance()->GetTime() << " ";

    // Loop over cells and find associated elements so in the same order as the cells in output files
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        unsigned elem_index = this->GetLocationIndexUsingCell(*cell_iter);

        // Hack that covers the case where the element is associated with a cell that has just been killed (#1129)
        bool elem_corresponds_to_dead_cell = false;

        if (this->IsCellAttachedToLocationIndex(elem_index))
        {
            elem_corresponds_to_dead_cell = this->GetCellUsingLocationIndex(elem_index)->IsDead();
        }

        // Write node data to file
        if (!(GetElement(elem_index)->IsDeleted()) && !elem_corresponds_to_dead_cell)
        {
           // Write element index to file
            *(this->mpCellVolumesFile) << elem_index << " ";

            // Write cell ID to file
            unsigned cell_index = cell_iter->GetCellId();
            *(this->mpCellVolumesFile) << cell_index << " ";

            // Write centroid location to file
            c_vector<double, DIM> centroid_location = mpPottsMesh->GetCentroidOfElement(elem_index);

            *(this->mpCellVolumesFile) << centroid_location[0] << " ";
            *(this->mpCellVolumesFile) << centroid_location[1] << " ";

            // Write cell volume (in 3D) or area (in 2D) to file
            double cell_volume = this->GetVolumeOfCell(*cell_iter);
            *(this->mpCellVolumesFile) << cell_volume << " ";
        }
    }
    *(this->mpCellVolumesFile) << "\n";
}

template<unsigned DIM>
void PottsBasedCellPopulation<DIM>::GenerateCellResultsAndWriteToFiles()
{
    // Set up cell type counter
    unsigned num_cell_types = this->mCellProliferativeTypeCount.size();
    std::vector<unsigned> cell_type_counter(num_cell_types);
    for (unsigned i=0; i<num_cell_types; i++)
    {
        cell_type_counter[i] = 0;
    }

    // Set up cell cycle phase counter
    unsigned num_cell_cycle_phases = this->mCellCyclePhaseCount.size();
    std::vector<unsigned> cell_cycle_phase_counter(num_cell_cycle_phases);
    for (unsigned i=0; i<num_cell_cycle_phases; i++)
    {
        cell_cycle_phase_counter[i] = 0;
    }

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        this->GenerateCellResults(*cell_iter, cell_type_counter, cell_cycle_phase_counter);
    }

    this->WriteCellResultsToFiles(cell_type_counter, cell_cycle_phase_counter);
}

template<unsigned DIM>
double PottsBasedCellPopulation<DIM>::GetWidth(const unsigned& rDimension)
{
    // Call GetWidth() on the mesh
    double width = this->mrMesh.GetWidth(rDimension);

    return width;
}

template<unsigned DIM>
void PottsBasedCellPopulation<DIM>::AddUpdateRule(boost::shared_ptr<AbstractPottsUpdateRule<DIM> > pUpdateRule)
{
    mUpdateRuleCollection.push_back(pUpdateRule);
}

template<unsigned DIM>
void PottsBasedCellPopulation<DIM>::RemoveAllUpdateRules()
{
	mUpdateRuleCollection.clear();
}

template<unsigned DIM>
const std::vector<boost::shared_ptr<AbstractPottsUpdateRule<DIM> > >& PottsBasedCellPopulation<DIM>::rGetUpdateRuleCollection() const
{
    return mUpdateRuleCollection;
}

template<unsigned DIM>
void PottsBasedCellPopulation<DIM>::CreateElementTessellation()
{
    ///\todo implement this method (#1666)
//  delete mpElementTessellation;
//
//    ///\todo this code would need to be extended if the domain were required to be periodic
//
//  std::vector<Node<2>*> nodes;
//  for (unsigned node_index=0; node_index<mrMesh.GetNumNodes(); node_index++)
//  {
//      Node<2>* p_temp_node = mrMesh.GetNode(node_index);
//      nodes.push_back(p_temp_node);
//  }
//  MutableMesh<2,2> mesh(nodes);
//    mpElementTessellation = new VertexMesh<2,2>(mesh);
}

template<unsigned DIM>
VertexMesh<DIM,DIM>* PottsBasedCellPopulation<DIM>::GetElementTessellation()
{
//    assert(mpElementTessellation != NULL);
    return mpElementTessellation;
}

template<unsigned DIM>
void PottsBasedCellPopulation<DIM>::CreateMutableMesh()
{
    delete mpMutableMesh;

    // Get the nodes of the PottsMesh
    std::vector<Node<DIM>*> nodes;
    for (unsigned node_index=0; node_index<this->mrMesh.GetNumNodes(); node_index++)
    {
      c_vector<double, DIM> location = this->mrMesh.GetNode(node_index)->rGetLocation();
      nodes.push_back(new Node<DIM>(node_index, location));
    }

    mpMutableMesh = new MutableMesh<DIM,DIM>(nodes);
}

template<unsigned DIM>
MutableMesh<DIM,DIM>* PottsBasedCellPopulation<DIM>::GetMutableMesh()
{
    assert(mpMutableMesh);
    return mpMutableMesh;
}

template<unsigned DIM>
void PottsBasedCellPopulation<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<Temperature>" << mTemperature << "</Temperature>\n";
    *rParamsFile << "\t\t<NumSweepsPerTimestep>" << mNumSweepsPerTimestep << "</NumSweepsPerTimestep>\n";

    // Call method on direct parent class
    AbstractOnLatticeCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}

template<unsigned DIM>
std::set<unsigned> PottsBasedCellPopulation<DIM>::GetNeighbouringNodeIndices(unsigned index)
{
    EXCEPTION("Cannot call GetNeighbouringNodeIndices() on a PottsBasedCellPopulation, need to go through the PottsMesh instead");
    std::set<unsigned> neighbouring_node_indices;
    return neighbouring_node_indices;
}

template<unsigned DIM>
void PottsBasedCellPopulation<DIM>::SetTemperature(double temperature)
{
    mTemperature = temperature;
}

template<unsigned DIM>
double PottsBasedCellPopulation<DIM>::GetTemperature()
{
    return mTemperature;
}

template<unsigned DIM>
void PottsBasedCellPopulation<DIM>::SetNumSweepsPerTimestep(unsigned numSweepsPerTimestep)
{
    mNumSweepsPerTimestep = numSweepsPerTimestep;
}

template<unsigned DIM>
unsigned PottsBasedCellPopulation<DIM>::GetNumSweepsPerTimestep()
{
    return mNumSweepsPerTimestep;
}

template<unsigned DIM>
void PottsBasedCellPopulation<DIM>::WriteVtkResultsToFile()
{
#ifdef CHASTE_VTK
    std::stringstream time;
    time << SimulationTime::Instance()->GetTimeStepsElapsed();
    VtkMeshWriter<DIM, DIM> mesh_writer(this->mDirPath, "results_"+time.str(), false);

    unsigned num_nodes = GetNumNodes();
    std::vector<double> cell_types;
    std::vector<double> cell_mutation_states;
    std::vector<double> cell_labels;
    std::vector<double> elem_ids;
    cell_types.reserve(num_nodes);
    cell_mutation_states.reserve(num_nodes);
    cell_labels.reserve(num_nodes);
    elem_ids.reserve(num_nodes);
    std::vector<std::vector<double> > cellwise_data;

    unsigned num_cell_data_items = 0;
    std::vector<std::string> cell_data_names;
    //We assume that the first cell is representative of all cells
    num_cell_data_items = this->Begin()->GetCellData()->GetNumItems();
    cell_data_names = this->Begin()->GetCellData()->GetKeys();

    for (unsigned var=0; var<num_cell_data_items; var++)
    {
        std::vector<double> cellwise_data_var(num_nodes);
        cellwise_data.push_back(cellwise_data_var);
    }
    for (typename AbstractMesh<DIM,DIM>::NodeIterator iter = mpPottsMesh->GetNodeIteratorBegin();
         iter != mpPottsMesh->GetNodeIteratorEnd();
         ++iter)
    {
        std::set<unsigned> element_indices = iter->rGetContainingElementIndices();

        if (element_indices.empty())
        {
            // No elements associated with this gridpoint
            cell_types.push_back(-1.0);
            elem_ids.push_back(-1.0);
            if (this->mOutputCellMutationStates)
            {
                cell_mutation_states.push_back(-1.0);
                cell_labels.push_back(-1.0);
            }
        }
        else
        {
            // The number of elements should be zero or one
            assert(element_indices.size() == 1);

            unsigned element_index = *(element_indices.begin());
            elem_ids.push_back((double)element_index);

            CellPtr p_cell = this->GetCellUsingLocationIndex(element_index);
            double cell_type = p_cell->GetCellCycleModel()->GetCellProliferativeType();
            cell_types.push_back(cell_type);

            if (this->mOutputCellMutationStates)
            {
                double cell_mutation_state = p_cell->GetMutationState()->GetColour();
                cell_mutation_states.push_back(cell_mutation_state);

                double cell_label = 0.0;
                if (p_cell->HasCellProperty<CellLabel>())
                {
                    CellPropertyCollection collection = p_cell->rGetCellPropertyCollection().GetProperties<CellLabel>();
                    boost::shared_ptr<CellLabel> p_label = boost::static_pointer_cast<CellLabel>(collection.GetProperty());
                    cell_label = p_label->GetColour();
                }
                cell_labels.push_back(cell_label);
            }
            for (unsigned var=0; var<num_cell_data_items; var++)
            {
                cellwise_data[var][iter->GetIndex()] = p_cell->GetCellData()->GetItem(cell_data_names[var]);
            }
        }
    }

    assert(cell_types.size() == num_nodes);
    assert(elem_ids.size() == num_nodes);

    mesh_writer.AddPointData("Element index", elem_ids);
    mesh_writer.AddPointData("Cell types", cell_types);

    if (this->mOutputCellMutationStates)
    {
        assert(cell_mutation_states.size() == num_nodes);
        mesh_writer.AddPointData("Mutation states", cell_mutation_states);
        assert(cell_labels.size() == num_nodes);
        mesh_writer.AddPointData("Cell labels", cell_labels);
    }

    if (num_cell_data_items > 0)
    {
        for (unsigned var=0; var<cellwise_data.size(); var++)
        {
            mesh_writer.AddPointData(cell_data_names[var], cellwise_data[var]);
        }
    }

    /*
     * The current VTK writer can only write things which inherit from AbstractTetrahedralMeshWriter.
     * For now, we do an explicit conversion to NodesOnlyMesh. This can be written to VTK then visualized as glyphs.
     */
    NodesOnlyMesh<DIM> temp_mesh;
    temp_mesh.ConstructNodesWithoutMesh(*mpPottsMesh);
    mesh_writer.WriteFilesUsingMesh(temp_mesh);

    *(this->mpVtkMetaFile) << "        <DataSet timestep=\"";
    *(this->mpVtkMetaFile) << time.str();
    *(this->mpVtkMetaFile) << "\" group=\"\" part=\"0\" file=\"results_";
    *(this->mpVtkMetaFile) << time.str();
    *(this->mpVtkMetaFile) << ".vtu\"/>\n";
#endif //CHASTE_VTK
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class PottsBasedCellPopulation<1>;
template class PottsBasedCellPopulation<2>;
template class PottsBasedCellPopulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PottsBasedCellPopulation)
