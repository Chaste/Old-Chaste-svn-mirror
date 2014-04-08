/*

Copyright (c) 2005-2014, University of Oxford.
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

#include "VertexBasedCellPopulation.hpp"
#include <boost/foreach.hpp>
#include "VertexMeshWriter.hpp"
#include "Warnings.hpp"
#include "ChasteSyscalls.hpp"
#include "IsNan.hpp"
#include "ShortAxisDivisionRule.hpp"

// Cell writers
#include "CellAgesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellVolumesWriter.hpp"

// Cell population writers
#include "CellMutationStatesCountWriter.hpp"
#include "CellPopulationElementWriter.hpp"
#include "VertexT1SwapLocationsWriter.hpp"
#include "VertexT2SwapLocationsWriter.hpp"
#include "VertexT3SwapLocationsWriter.hpp"

template<unsigned DIM>
VertexBasedCellPopulation<DIM>::VertexBasedCellPopulation(MutableVertexMesh<DIM, DIM>& rMesh,
                                          std::vector<CellPtr>& rCells,
                                          bool deleteMesh,
                                          bool validate,
                                          const std::vector<unsigned> locationIndices)
    : AbstractOffLatticeCellPopulation<DIM>(rMesh, rCells, locationIndices),
      mDeleteMesh(deleteMesh),
      mOutputCellRearrangementLocations(true)
{
    mpMutableVertexMesh = static_cast<MutableVertexMesh<DIM, DIM>* >(&(this->mrMesh));
    mpDivisionRule.reset(new ShortAxisDivisionRule<DIM>());

    // If no location indices are specified, associate with elements from the mesh (assumed to be sequentially ordered).
    std::list<CellPtr>::iterator it = this->mCells.begin();
    for (unsigned i=0; it != this->mCells.end(); ++it, ++i)
    {
        unsigned index = locationIndices.empty() ? i : locationIndices[i]; // assume that the ordering matches
        AbstractCellPopulation<DIM, DIM>::AddCellUsingLocationIndex(index,*it);
    }

    // Check each element has only one cell attached
    if (validate)
    {
        Validate();
    }
}

template<unsigned DIM>
VertexBasedCellPopulation<DIM>::VertexBasedCellPopulation(MutableVertexMesh<DIM, DIM>& rMesh)
    : AbstractOffLatticeCellPopulation<DIM>(rMesh),
      mDeleteMesh(true),
      mOutputCellRearrangementLocations(true)
{
    mpMutableVertexMesh = static_cast<MutableVertexMesh<DIM, DIM>* >(&(this->mrMesh));
}

template<unsigned DIM>
VertexBasedCellPopulation<DIM>::~VertexBasedCellPopulation()
{
    if (mDeleteMesh)
    {
        delete &this->mrMesh;
    }
}

template<unsigned DIM>
double VertexBasedCellPopulation<DIM>::GetDampingConstant(unsigned nodeIndex)
{
    // Take the average of the cells containing this vertex
    double average_damping_constant = 0.0;

    std::set<unsigned> containing_elements = GetNode(nodeIndex)->rGetContainingElementIndices();
    double temp = 1.0/((double) containing_elements.size());

    for (std::set<unsigned>::iterator iter = containing_elements.begin();
         iter != containing_elements.end();
         ++iter)
    {
        CellPtr p_cell = this->GetCellUsingLocationIndex(*iter);
        bool cell_is_wild_type = p_cell->GetMutationState()->IsType<WildTypeCellMutationState>();
        bool cell_is_labelled = p_cell->HasCellProperty<CellLabel>();

        if (cell_is_wild_type && !cell_is_labelled)
        {
            average_damping_constant += this->GetDampingConstantNormal()*temp;
        }
        else
        {
            average_damping_constant += this->GetDampingConstantMutant()*temp;
        }
    }

    return average_damping_constant;
}

template<unsigned DIM>
MutableVertexMesh<DIM, DIM>& VertexBasedCellPopulation<DIM>::rGetMesh()
{
    return *mpMutableVertexMesh;
}

template<unsigned DIM>
const MutableVertexMesh<DIM, DIM>& VertexBasedCellPopulation<DIM>::rGetMesh() const
{
    return *mpMutableVertexMesh;
}

template<unsigned DIM>
VertexElement<DIM, DIM>* VertexBasedCellPopulation<DIM>::GetElement(unsigned elementIndex)
{
    return mpMutableVertexMesh->GetElement(elementIndex);
}

template<unsigned DIM>
unsigned VertexBasedCellPopulation<DIM>::GetNumNodes()
{
    return this->mrMesh.GetNumNodes();
}

template<unsigned DIM>
c_vector<double, DIM> VertexBasedCellPopulation<DIM>::GetLocationOfCellCentre(CellPtr pCell)
{
    return mpMutableVertexMesh->GetCentroidOfElement(this->mCellLocationMap[pCell.get()]);
}

template<unsigned DIM>
Node<DIM>* VertexBasedCellPopulation<DIM>::GetNode(unsigned index)
{
    return this->mrMesh.GetNode(index);
}

template<unsigned DIM>
unsigned VertexBasedCellPopulation<DIM>::AddNode(Node<DIM>* pNewNode)
{
    return mpMutableVertexMesh->AddNode(pNewNode);
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::SetNode(unsigned nodeIndex, ChastePoint<DIM>& rNewLocation)
{
    mpMutableVertexMesh->SetNode(nodeIndex, rNewLocation);
}

template<unsigned DIM>
VertexElement<DIM, DIM>* VertexBasedCellPopulation<DIM>::GetElementCorrespondingToCell(CellPtr pCell)
{
    return mpMutableVertexMesh->GetElement(this->GetLocationIndexUsingCell(pCell));
}

template<unsigned DIM>
unsigned VertexBasedCellPopulation<DIM>::GetNumElements()
{
    return mpMutableVertexMesh->GetNumElements();
}

template<unsigned DIM>
CellPtr VertexBasedCellPopulation<DIM>::AddCell(CellPtr pNewCell,
                                                const c_vector<double,DIM>& rCellDivisionVector,
                                                CellPtr pParentCell)
{
    // Get the element associated with this cell
    VertexElement<DIM, DIM>* p_element = GetElementCorrespondingToCell(pParentCell);

    // Divide the element
    unsigned new_element_index = mpMutableVertexMesh->DivideElementAlongGivenAxis(p_element,
                                                                                  rCellDivisionVector,
                                                                                  true);
    // Associate the new cell with the element
    this->mCells.push_back(pNewCell);

    // Update location cell map
    CellPtr p_created_cell = this->mCells.back();
    this->SetCellUsingLocationIndex(new_element_index,p_created_cell);
    this->mCellLocationMap[p_created_cell.get()] = new_element_index;
    return p_created_cell;
}

template<unsigned DIM>
unsigned VertexBasedCellPopulation<DIM>::RemoveDeadCells()
{
    unsigned num_removed = 0;

    for (std::list<CellPtr>::iterator it = this->mCells.begin();
         it != this->mCells.end();
         )
    {
        if ((*it)->IsDead())
        {
            // Count the cell as dead
            num_removed++;

            // Remove the element from the mesh if it is not deleted yet
            ///\todo (#2489) this should cause an error - we should fix this!
            if (!(this->GetElement(this->GetLocationIndexUsingCell((*it)))->IsDeleted()))
            {
                // This warning relies on the fact that there is only one other possibility for
                // vertex elements to be marked as deleted: a T2 swap
                WARNING("A Cell is removed without performing a T2 swap. This leaves a void in the mesh.");
                mpMutableVertexMesh->DeleteElementPriorToReMesh(this->GetLocationIndexUsingCell((*it)));
            }

            // Delete the cell
            it = this->mCells.erase(it);
        }
        else
        {
            ++it;
        }
    }
    return num_removed;
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::UpdateNodeLocations(double dt)
{
    // Iterate over all nodes associated with real cells to update their positions
    for (unsigned node_index=0; node_index<GetNumNodes(); node_index++)
    {
        // Get the damping constant for this node
        double damping_const = this->GetDampingConstant(node_index);

        // Compute the displacement of this node
        c_vector<double, DIM> displacement = dt*this->GetNode(node_index)->rGetAppliedForce()/damping_const;

        /*
         * If the displacement of this node is greater than half the cell rearrangement threshold,
         * this could result in nodes moving into the interior of other elements, which should not
         * be possible. Therefore in this case we restrict the displacement of the node to the cell
         * rearrangement threshold and warn the user that a smaller timestep should be used. This
         * restriction ensures that vertex elements remain well defined (see #1376).
         */
        if (norm_2(displacement) > 0.5*mpMutableVertexMesh->GetCellRearrangementThreshold())
        {
            WARN_ONCE_ONLY("Vertices are moving more than half the CellRearrangementThreshold. This could cause elements to become inverted so the motion has been restricted. Use a smaller timestep to avoid these warnings.");
            displacement *= 0.5*mpMutableVertexMesh->GetCellRearrangementThreshold()/norm_2(displacement);
        }

        // Get new node location
        c_vector<double, DIM> new_node_location = this->GetNode(node_index)->rGetLocation() + displacement;

        for (unsigned i=0; i<DIM; i++)
        {
            assert(!std::isnan(new_node_location(i)));
        }

        // Create ChastePoint for new node location
        ChastePoint<DIM> new_point(new_node_location);

        // Move the node
        this->SetNode(node_index, new_point);
    }
}

template<unsigned DIM>
bool VertexBasedCellPopulation<DIM>::IsCellAssociatedWithADeletedLocation(CellPtr pCell)
{
    return GetElementCorrespondingToCell(pCell)->IsDeleted();;
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::Update(bool hasHadBirthsOrDeaths)
{
    VertexElementMap element_map(mpMutableVertexMesh->GetNumAllElements());

    mpMutableVertexMesh->ReMesh(element_map);

    if (!element_map.IsIdentityMap())
    {
        // Fix up the mappings between CellPtrs and VertexElements
        ///\todo We want to make these maps private, so we need a better way of doing the code below.
        std::map<Cell*, unsigned> old_map = this->mCellLocationMap;

        this->mCellLocationMap.clear();
        this->mLocationCellMap.clear();

        for (std::list<CellPtr>::iterator cell_iter = this->mCells.begin();
             cell_iter != this->mCells.end();
             ++cell_iter)
        {
            // The cell vector should only ever contain living cells
            unsigned old_elem_index = old_map[(*cell_iter).get()];
            assert(!element_map.IsDeleted(old_elem_index));

            unsigned new_elem_index = element_map.GetNewIndex(old_elem_index);
            this->SetCellUsingLocationIndex(new_elem_index, *cell_iter);
        }

        // Check that each VertexElement has only one CellPtr associated with it in the updated cell population
        Validate();
    }

    element_map.ResetToIdentity();
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::Validate()
{
    // Check each element has only one cell attached
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
            // This should never be reached as you can only set one cell per element index
            EXCEPTION("Element " << i << " appears to have " << validated_element[i] << " cells associated with it");
        }
    }
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::AcceptPopulationWriter(boost::shared_ptr<AbstractCellPopulationWriter<DIM, DIM> > pPopulationWriter)
{
    pPopulationWriter->Visit(this);
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::AcceptCellWriter(boost::shared_ptr<AbstractCellWriter<DIM, DIM> > pCellWriter, CellPtr pCell)
{
    pCellWriter->VisitCell(pCell, this);
}

template<unsigned DIM>
double VertexBasedCellPopulation<DIM>::GetVolumeOfCell(CellPtr pCell)
{
    // Get the vertex element index corresponding to this cell
    unsigned elem_index = this->GetLocationIndexUsingCell(pCell);

    // Get the cell's volume from the vertex mesh
    double cell_volume = mpMutableVertexMesh->GetVolumeOfElement(elem_index);

    return cell_volume;
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::WriteVtkResultsToFile(const std::string& rDirectory)
{
#ifdef CHASTE_VTK
    SimulationTime* p_time = SimulationTime::Instance();

    VertexMeshWriter<DIM, DIM> mesh_writer(rDirectory, "results", false);
    std::stringstream time;
    time << p_time->GetTimeStepsElapsed();

    unsigned num_elements = mpMutableVertexMesh->GetNumElements();
    std::vector<double> cell_types(num_elements);
    std::vector<double> cell_labels(num_elements);
    std::vector<double> cell_ancestors(num_elements);
    std::vector<double> cell_mutation_states(num_elements);
    std::vector<double> cell_ages(num_elements);
    std::vector<double> cell_cycle_phases(num_elements);
    std::vector<double> cell_volumes(num_elements);
    std::vector<std::vector<double> > cellwise_data;

    unsigned num_cell_data_items = 0;
    std::vector<std::string> cell_data_names;

    // We assume that the first cell is representative of all cells
    num_cell_data_items = this->Begin()->GetCellData()->GetNumItems();
    cell_data_names = this->Begin()->GetCellData()->GetKeys();

    for (unsigned var=0; var<num_cell_data_items; var++)
    {
        std::vector<double> cellwise_data_var(num_elements);
        cellwise_data.push_back(cellwise_data_var);
    }

    // Loop over vertex elements
    for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = mpMutableVertexMesh->GetElementIteratorBegin();
         elem_iter != mpMutableVertexMesh->GetElementIteratorEnd();
         ++elem_iter)
    {
        // Get index of this element in the vertex mesh
        unsigned elem_index = elem_iter->GetIndex();

        // Get the cell corresponding to this element
        CellPtr p_cell = this->GetCellUsingLocationIndex(elem_index);
        assert(p_cell);

        double cell_label = 0.0;
        if (p_cell->HasCellProperty<CellLabel>())
        {
            CellPropertyCollection collection = p_cell->rGetCellPropertyCollection().GetProperties<CellLabel>();
            boost::shared_ptr<CellLabel> p_label = boost::static_pointer_cast<CellLabel>(collection.GetProperty());
            cell_label = p_label->GetColour();
        }
        cell_labels[elem_index] = cell_label;

        if (this-> template HasWriter<CellAncestorWriter>())
        {
            double ancestor_index = (p_cell->GetAncestor() == UNSIGNED_UNSET) ? (-1.0) : (double)p_cell->GetAncestor();
            cell_ancestors[elem_index] = ancestor_index;
        }
        if (this-> template HasWriter<CellProliferativeTypesWriter>())
        {
            double cell_type = p_cell->GetCellProliferativeType()->GetColour();
            cell_types[elem_index] = cell_type;
        }
        if (this-> template HasWriter<CellMutationStatesCountWriter>())
        {
            double mutation_state = p_cell->GetMutationState()->GetColour();
            cell_mutation_states[elem_index] = mutation_state;
        }
        if (this-> template HasWriter<CellAgesWriter>())
        {
            double age = p_cell->GetAge();
            cell_ages[elem_index] = age;
        }
        if (this-> template HasWriter<CellProliferativePhasesWriter>())
        {
            double cycle_phase = p_cell->GetCellCycleModel()->GetCurrentCellCyclePhase();
            cell_cycle_phases[elem_index] = cycle_phase;
        }
        if (this-> template HasWriter<CellVolumesWriter>())
        {
            double cell_volume = mpMutableVertexMesh->GetVolumeOfElement(elem_index);
            cell_volumes[elem_index] = cell_volume;
        }
        for (unsigned var=0; var<num_cell_data_items; var++)
        {
            cellwise_data[var][elem_index] = p_cell->GetCellData()->GetItem(cell_data_names[var]);
        }
    }

    mesh_writer.AddCellData("Cell labels", cell_labels);
    if (this-> template HasWriter<CellProliferativeTypesWriter>())
    {
        mesh_writer.AddCellData("Cell types", cell_types);
    }
    if (this-> template HasWriter<CellAncestorWriter>())
    {
        mesh_writer.AddCellData("Ancestors", cell_ancestors);
    }
    if (this-> template HasWriter<CellMutationStatesCountWriter>())
    {
        mesh_writer.AddCellData("Mutation states", cell_mutation_states);
    }
    if (this-> template HasWriter<CellAgesWriter>())
    {
        mesh_writer.AddCellData("Ages", cell_ages);
    }
    if (this-> template HasWriter<CellProliferativePhasesWriter>())
    {
        mesh_writer.AddCellData("Cycle phases", cell_cycle_phases);
    }
    if (this-> template HasWriter<CellVolumesWriter>())
    {
        mesh_writer.AddCellData("Cell volumes", cell_volumes);
    }
    if (num_cell_data_items > 0)
    {
        for (unsigned var=0; var<cellwise_data.size(); var++)
        {
            mesh_writer.AddCellData(cell_data_names[var], cellwise_data[var]);
        }
    }

    mesh_writer.WriteVtkUsingMesh(*mpMutableVertexMesh, time.str());
    *(this->mpVtkMetaFile) << "        <DataSet timestep=\"";
    *(this->mpVtkMetaFile) << p_time->GetTimeStepsElapsed();
    *(this->mpVtkMetaFile) << "\" group=\"\" part=\"0\" file=\"results_";
    *(this->mpVtkMetaFile) << p_time->GetTimeStepsElapsed();
    *(this->mpVtkMetaFile) << ".vtu\"/>\n";
#endif //CHASTE_VTK
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::OpenWritersFiles(const std::string& rDirectory)
{
    if (this->mOutputResultsForChasteVisualizer)
    {
        if (!this-> template HasWriter<CellPopulationElementWriter>())
        {
            this-> template AddPopulationWriter<CellPopulationElementWriter>();
        }
    }

    if (mOutputCellRearrangementLocations)
    {
        if (!this-> template HasWriter<VertexT1SwapLocationsWriter>())
        {
            this-> template AddPopulationWriter<VertexT1SwapLocationsWriter>();
        }

        if (!this-> template HasWriter<VertexT2SwapLocationsWriter>())
        {
            this-> template AddPopulationWriter<VertexT2SwapLocationsWriter>();
        }

        if (!this-> template HasWriter<VertexT3SwapLocationsWriter>())
        {
            this-> template AddPopulationWriter<VertexT3SwapLocationsWriter>();
        }
    }

    AbstractCellPopulation<DIM>::OpenWritersFiles(rDirectory);
}

template<unsigned DIM>
bool VertexBasedCellPopulation<DIM>::GetOutputCellRearrangementLocations()
{
    return mOutputCellRearrangementLocations;
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::SetOutputCellRearrangementLocations(bool outputCellRearrangementLocations)
{
    mOutputCellRearrangementLocations = outputCellRearrangementLocations;
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<CellRearrangementThreshold>" << mpMutableVertexMesh->GetCellRearrangementThreshold() << "</CellRearrangementThreshold>\n";
    *rParamsFile << "\t\t<T2Threshold>" <<  mpMutableVertexMesh->GetT2Threshold() << "</T2Threshold>\n";
    *rParamsFile << "\t\t<CellRearrangementRatio>" << mpMutableVertexMesh->GetCellRearrangementRatio() << "</CellRearrangementRatio>\n";
    *rParamsFile << "\t\t<OutputCellRearrangementLocations>" << mOutputCellRearrangementLocations << "</OutputCellRearrangementLocations>\n";

    //Add the division rule parameters
    *rParamsFile << "\t\t<DivisionRule>\n";
    mpDivisionRule->OutputCellDivisionRuleInfo(rParamsFile);
    *rParamsFile << "\t\t</DivisionRule>\n";


    // Call method on direct parent class
    AbstractOffLatticeCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}

template<unsigned DIM>
double VertexBasedCellPopulation<DIM>::GetWidth(const unsigned& rDimension)
{
    // Call GetWidth() on the mesh
    double width = this->mrMesh.GetWidth(rDimension);

    return width;
}

template<unsigned DIM>
std::set<unsigned> VertexBasedCellPopulation<DIM>::GetNeighbouringNodeIndices(unsigned index)
{
    return mpMutableVertexMesh->GetNeighbouringNodeIndices(index);
}

template<unsigned DIM>
boost::shared_ptr<AbstractCellDivisionRule<DIM> > VertexBasedCellPopulation<DIM>::GetDivisionRule()
{
    return mpDivisionRule;
}

template<unsigned DIM>
void VertexBasedCellPopulation<DIM>::SetDivisionRule(boost::shared_ptr<AbstractCellDivisionRule<DIM> > pDivisionRule)
{
    mpDivisionRule = pDivisionRule;
}

template<unsigned DIM>
TetrahedralMesh<DIM, DIM>* VertexBasedCellPopulation<DIM>::GetTetrahedralMeshUsingVertexMesh()
{
    // This method only works in 2D sequential
    assert(DIM == 2);
    assert(PetscTools::IsSequential());

    unsigned num_vertex_nodes = mpMutableVertexMesh->GetNumNodes();
    unsigned num_vertex_elements = mpMutableVertexMesh->GetNumElements();

    std::string mesh_file_name = "mesh";

    // Get a unique temporary foldername
    std::stringstream pid;
    pid << getpid();
    OutputFileHandler output_file_handler("2D_temporary_tetrahedral_mesh_" + pid.str());
    std::string output_dir = output_file_handler.GetOutputDirectoryFullPath();

    // Compute the number of nodes in the TetrahedralMesh
    unsigned num_tetrahedral_nodes = num_vertex_nodes + num_vertex_elements;

    // Write node file
    out_stream p_node_file = output_file_handler.OpenOutputFile(mesh_file_name+".node");
    (*p_node_file) << std::scientific;
    (*p_node_file) << std::setprecision(20);
    (*p_node_file) << num_tetrahedral_nodes << "\t2\t0\t1" << std::endl;

    // Begin by writing each node in the VertexMesh
    for (unsigned node_index=0; node_index<num_vertex_nodes; node_index++)
    {
        Node<DIM>* p_node = mpMutableVertexMesh->GetNode(node_index);

        ///\todo will the nodes in mpMutableVertexMesh always have indices 0,1,2,...? (#2221)
        unsigned index = p_node->GetIndex();

        c_vector<double, DIM> location = p_node->rGetLocation();

        ///\todo we do not yet consider boundaryness in vertex meshes (see #1558)
        unsigned is_boundary_node = p_node->IsBoundaryNode() ? 1 : 0;

        (*p_node_file) << index << "\t" << location[0] << "\t" << location[1] << "\t" << is_boundary_node << std::endl;
    }

    // Now write an additional node at each VertexElement's centroid
    unsigned num_tetrahedral_elements = 0;
    for (unsigned vertex_elem_index=0; vertex_elem_index<num_vertex_elements; vertex_elem_index++)
    {
        unsigned index = num_vertex_nodes + vertex_elem_index;

        c_vector<double, DIM> location = mpMutableVertexMesh->GetCentroidOfElement(vertex_elem_index);

        // Any node located at a VertexElement's centroid will not be a boundary node
        unsigned is_boundary_node = 0;
        (*p_node_file) << index << "\t" << location[0] << "\t" << location[1] << "\t" << is_boundary_node << std::endl;

        // Also keep track of how many tetrahedral elements there will be
        num_tetrahedral_elements += mpMutableVertexMesh->GetElement(vertex_elem_index)->GetNumNodes();
    }
    p_node_file->close();

    // Write element file
    out_stream p_elem_file = output_file_handler.OpenOutputFile(mesh_file_name+".ele");
    (*p_elem_file) << std::scientific;
    (*p_elem_file) << num_tetrahedral_elements << "\t3\t0" << std::endl;

    std::set<std::pair<unsigned, unsigned> > tetrahedral_edges;

    unsigned tetrahedral_elem_index = 0;
    for (unsigned vertex_elem_index=0; vertex_elem_index<num_vertex_elements; vertex_elem_index++)
    {
        VertexElement<DIM, DIM>* p_vertex_element = mpMutableVertexMesh->GetElement(vertex_elem_index);

        // Iterate over nodes owned by this VertexElement
        unsigned num_nodes_in_vertex_element = p_vertex_element->GetNumNodes();
        for (unsigned local_index=0; local_index<num_nodes_in_vertex_element; local_index++)
        {
            unsigned node_0_index = p_vertex_element->GetNodeGlobalIndex(local_index);
            unsigned node_1_index = p_vertex_element->GetNodeGlobalIndex((local_index+1)%num_nodes_in_vertex_element);
            unsigned node_2_index = num_vertex_nodes + vertex_elem_index;

            (*p_elem_file) << tetrahedral_elem_index++ << "\t" << node_0_index << "\t" << node_1_index << "\t" << node_2_index << std::endl;

            // Add edges to the set if they are not already present
            std::pair<unsigned, unsigned> edge_0 = this->CreateOrderedPair(node_0_index, node_1_index);
            std::pair<unsigned, unsigned> edge_1 = this->CreateOrderedPair(node_1_index, node_2_index);
            std::pair<unsigned, unsigned> edge_2 = this->CreateOrderedPair(node_2_index, node_0_index);

            tetrahedral_edges.insert(edge_0);
            tetrahedral_edges.insert(edge_1);
            tetrahedral_edges.insert(edge_2);
        }
    }
    p_elem_file->close();

    // Write edge file
    out_stream p_edge_file = output_file_handler.OpenOutputFile(mesh_file_name+".edge");
    (*p_edge_file) << std::scientific;
    (*p_edge_file) << tetrahedral_edges.size() << "\t1" << std::endl;

    unsigned edge_index = 0;
    for (std::set<std::pair<unsigned, unsigned> >::iterator edge_iter = tetrahedral_edges.begin();
         edge_iter != tetrahedral_edges.end();
         ++edge_iter)
    {
        std::pair<unsigned, unsigned> this_edge;

        ///\todo we do not yet consider boundaryness in vertex meshes (see #1558)
        (*p_edge_file) << edge_index++ << "\t" << this_edge.first << "\t" << this_edge.second << "\t" << 0 << std::endl;
    }
    p_edge_file->close();

    // Having written the mesh to file, now construct it using TrianglesMeshReader
    TetrahedralMesh<DIM, DIM>* p_mesh = new TetrahedralMesh<DIM, DIM>;
    // Nested scope so reader is destroyed before we remove the temporary files.
    {
        TrianglesMeshReader<DIM, DIM> mesh_reader(output_dir + mesh_file_name);
        p_mesh->ConstructFromMeshReader(mesh_reader);
    }

    // Delete the temporary files
    output_file_handler.FindFile("").Remove();

    // The original files have been deleted, it is better if the mesh object forgets about them
    p_mesh->SetMeshHasChangedSinceLoading();

    return p_mesh;
}

// Explicit instantiation
template class VertexBasedCellPopulation<1>;
template class VertexBasedCellPopulation<2>;
template class VertexBasedCellPopulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VertexBasedCellPopulation)
