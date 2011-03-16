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
#include "NodeBasedCellPopulation.hpp"
#include "CellwiseData.hpp"
#include "VtkMeshWriter.hpp"

template<unsigned DIM>
NodeBasedCellPopulation<DIM>::NodeBasedCellPopulation(const std::vector<Node<DIM>* > nodes,
                                      std::vector<CellPtr>& rCells,
                                      const std::vector<unsigned> locationIndices,
                                      bool deleteNodes)
    : AbstractCentreBasedCellPopulation<DIM>(rCells, locationIndices),
      mNodes(nodes.begin(), nodes.end()),
      mAddedNodes(true),
      mpBoxCollection(NULL),
      mDeleteNodes(deleteNodes),
      mMechanicsCutOffLength(DBL_MAX)
{
    Validate();
}

// archiving constructor
template<unsigned DIM>
NodeBasedCellPopulation<DIM>::NodeBasedCellPopulation(const std::vector<Node<DIM>* > nodes, double mechanicsCutOffLength, bool deleteNodes)
    : AbstractCentreBasedCellPopulation<DIM>(),
      mNodes(nodes.begin(), nodes.end()),
      mAddedNodes(true),
      mpBoxCollection(NULL),
      mDeleteNodes(deleteNodes),
      mMechanicsCutOffLength(mechanicsCutOffLength)
{
    // No Validate() because the cells are not associated with the cell population yet in archiving
}

template<unsigned DIM>
NodeBasedCellPopulation<DIM>::NodeBasedCellPopulation(const AbstractMesh<DIM,DIM>& rMesh,
                                      std::vector<CellPtr>& rCells)
    : AbstractCentreBasedCellPopulation<DIM>(rCells),
      mAddedNodes(false),
      mpBoxCollection(NULL),
      mDeleteNodes(true),
      mMechanicsCutOffLength(DBL_MAX)
{
    unsigned num_nodes = rMesh.GetNumNodes();

    mNodes.reserve(num_nodes);
    // Copy the actual node objects from mesh to this (mesh-less) cell population.
    for (unsigned i=0; i<num_nodes; i++)
    {
        Node<DIM>* p_node = new Node<DIM>(*(rMesh.GetNode(i)));
        mNodes.push_back(p_node);
    }
    mAddedNodes = true;
    Validate();
}

template<unsigned DIM>
NodeBasedCellPopulation<DIM>::~NodeBasedCellPopulation()
{
    Clear();

    // Free node memory
    if (mDeleteNodes)
    {
        for (unsigned i=0; i<mNodes.size(); i++)
        {
            delete mNodes[i];
        }
    }
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::Clear()
{
    delete mpBoxCollection;
    mpBoxCollection = NULL;
    mNodePairs.clear();
    mDeletedNodeIndices.clear();
    mAddedNodes = false;
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::Validate()
{
    std::vector<bool> validated_node(GetNumNodes());
    for (unsigned i=0; i<validated_node.size(); i++)
    {
        validated_node[i] = false;
    }

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        unsigned node_index = this->mCellLocationMap[(*cell_iter).get()];
        validated_node[node_index] = true;
    }

    for (unsigned i=0; i<validated_node.size(); i++)
    {
        if (!validated_node[i])
        {
            std::stringstream ss;
            ss << "Node " << i << " does not appear to have a cell associated with it";
            EXCEPTION(ss.str());
        }
    }
}

template<unsigned DIM>
std::vector<Node<DIM>* >& NodeBasedCellPopulation<DIM>::rGetNodes()
{
    return mNodes;
}

template<unsigned DIM>
const std::vector<Node<DIM>* >& NodeBasedCellPopulation<DIM>::rGetNodes() const
{
    return mNodes;
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::SplitUpIntoBoxes(double cutOffLength, c_vector<double, 2*DIM> domainSize)
{
    mpBoxCollection = new BoxCollection<DIM>(cutOffLength, domainSize);
    mpBoxCollection->SetupLocalBoxesHalfOnly();

    for (unsigned i=0; i<mNodes.size(); i++)
    {
        unsigned box_index = mpBoxCollection->CalculateContainingBox(mNodes[i]);
        mpBoxCollection->rGetBox(box_index).AddNode(mNodes[i]);
    }
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::FindMaxAndMin()
{
    c_vector<double, DIM> min_posn;
    c_vector<double, DIM> max_posn;
    for (unsigned i=0; i<DIM; i++)
    {
        min_posn(i) = DBL_MAX;
        max_posn(i) = -DBL_MAX;
    }

    for (unsigned i=0; i<mNodes.size(); i++)
    {
        c_vector<double, DIM> posn = this->GetNode(i)->rGetLocation();

        for (unsigned j=0; j<DIM; j++)
        {
            if (posn(j) > max_posn(j))
            {
                max_posn(j) = posn(j);
            }
            if (posn(j) < min_posn(j))
            {
                min_posn(j) = posn(j);
            }
        }
    }

    for (unsigned i=0; i<DIM; i++)
    {
        assert(min_posn(i) != DBL_MAX);
        mMinSpatialPositions(i) = min_posn(i);

        assert(max_posn(i) != -DBL_MAX);
        mMaxSpatialPositions(i) = max_posn(i);
    }
}

template<unsigned DIM>
Node<DIM>* NodeBasedCellPopulation<DIM>::GetNode(unsigned index)
{
    return mNodes[index];
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::SetNode(unsigned nodeIndex, ChastePoint<DIM>& rNewLocation)
{
    mNodes[nodeIndex]->SetPoint(rNewLocation);
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::Update(bool hasHadBirthsOrDeaths)
{
    if (hasHadBirthsOrDeaths)
    {
        // Create and reserve space for a temporary vector
        std::vector<Node<DIM>*> old_nodes;
        old_nodes.reserve(mNodes.size());

        // Store all non-deleted nodes in the temporary vector
        for (unsigned i=0; i<mNodes.size(); i++)
        {
            if (!mNodes[i]->IsDeleted())
            {
                old_nodes.push_back(mNodes[i]);
            }
            else
            {
                // Free node memory
                delete mNodes[i];
            }
        }

        std::map<unsigned,CellPtr> old_map = this->mLocationCellMap;
        mNodes.clear();

        // Clear maps
        this->mLocationCellMap.clear();
        this->mCellLocationMap.clear();

        // Update mNodes to new indices which go from 0 to NumNodes-1
        for (unsigned i=0; i<old_nodes.size(); i++)
        {
            // Get the living cell associated with the old node
            CellPtr p_live_cell = old_map[old_nodes[i]->GetIndex()];

            // Set the node up
            mNodes.push_back(old_nodes[i]);
            mNodes[i]->SetIndex(i);

            // Set the maps up
            this->mLocationCellMap[i] = p_live_cell;
            this->mCellLocationMap[p_live_cell.get()] = i;
        }

        // Remove current dead indices data
        Clear();

        Validate();
    }

    if (mpBoxCollection!=NULL)
    {
        delete mpBoxCollection;
    }

    FindMaxAndMin();

    // Something here to set up the domain size (max and min of each node position dimension)
    c_vector<double, 2*DIM> domain_size;

    for (unsigned i=0; i<DIM; i++)
    {
        domain_size(2*i) = mMinSpatialPositions(i);
        domain_size(2*i+1) = mMaxSpatialPositions(i);
    }

    if (mMechanicsCutOffLength == DBL_MAX)
    {
        std::string error =  std::string("NodeBasedCellPopulation cannot create boxes if the cut-off length has not been set - ")
                           + std::string("Call SetMechanicsCutOffLength on the CellPopulation ensuring it is larger than GetCutOffLength() on the force law");
        EXCEPTION(error);
    }

    // Add this parameter and suggest that mechanics systems set it.
    // Allocates memory for mpBoxCollection and does the splitting and putting nodes into boxes
    SplitUpIntoBoxes(mMechanicsCutOffLength, domain_size);

    mpBoxCollection->CalculateNodePairs(mNodes, mNodePairs);
}

template<unsigned DIM>
unsigned NodeBasedCellPopulation<DIM>::RemoveDeadCells()
{
    unsigned num_removed = 0;

    for (std::list<CellPtr>::iterator cell_iter = this->mCells.begin();
         cell_iter != this->mCells.end();
         ++cell_iter)
    {
        if ((*cell_iter)->IsDead())
        {
            // Remove the node
            num_removed++;
            this->GetNodeCorrespondingToCell(*cell_iter)->MarkAsDeleted();
            mDeletedNodeIndices.push_back(this->mCellLocationMap[(*cell_iter).get()]);
            cell_iter = this->mCells.erase(cell_iter);
            --cell_iter;
        }
    }
    return num_removed;
}

template<unsigned DIM>
unsigned NodeBasedCellPopulation<DIM>::AddNode(Node<DIM>* pNewNode)
{
    if (mDeletedNodeIndices.empty())
    {
        pNewNode->SetIndex(mNodes.size());
        mNodes.push_back(pNewNode);
    }
    else
    {
        unsigned index = mDeletedNodeIndices.back();
        pNewNode->SetIndex(index);
        mDeletedNodeIndices.pop_back();
        delete mNodes[index];
        mNodes[index] = pNewNode;
    }
    mAddedNodes = true;
    return pNewNode->GetIndex();
}

template<unsigned DIM>
unsigned NodeBasedCellPopulation<DIM>::GetNumNodes()
{
    return mNodes.size() - mDeletedNodeIndices.size();
}

template<unsigned DIM>
BoxCollection<DIM>* NodeBasedCellPopulation<DIM>::GetBoxCollection()
{
    return mpBoxCollection;
}

template<unsigned DIM>
std::set< std::pair<Node<DIM>*, Node<DIM>* > >& NodeBasedCellPopulation<DIM>::rGetNodePairs()
{
    if (mNodePairs.empty())
    {
        EXCEPTION("No node pairs set up, rGetNodePairs probably called before Update");
    }
    return mNodePairs;
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<MechanicsCutOffLength>"<< mMechanicsCutOffLength << "</MechanicsCutOffLength>\n";

    // Currently no specific parameters to output all come from parent classes

	// Call direct parent class method
	AbstractCentreBasedCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::SetMechanicsCutOffLength(double mechanicsCutOffLength)
{
    assert(mechanicsCutOffLength > 0.0);
    mMechanicsCutOffLength = mechanicsCutOffLength;
}

template<unsigned DIM>
double NodeBasedCellPopulation<DIM>::GetMechanicsCutOffLength()
{
    return mMechanicsCutOffLength;
}

template<unsigned DIM>
double NodeBasedCellPopulation<DIM>::GetWidth(const unsigned& rDimension)
{
    // Update the member variables mMinSpatialPositions and mMaxSpatialPositions
    FindMaxAndMin();

    // Compute the maximum distance between any nodes in this dimension
    double width = mMaxSpatialPositions(rDimension) - mMinSpatialPositions(rDimension);

    return width;
}

template<unsigned DIM>
void NodeBasedCellPopulation<DIM>::WriteVtkResultsToFile()
{
#ifdef CHASTE_VTK
    
    std::stringstream time;
    time << SimulationTime::Instance()->GetTimeStepsElapsed();
    VtkMeshWriter<DIM, DIM> mesh_writer(this->mDirPath, "results_"+time.str(), false);

    unsigned num_nodes = GetNumNodes();
    std::vector<double> cell_types(num_nodes);
    std::vector<double> cell_ancestors(num_nodes);
    std::vector<double> cell_mutation_states(num_nodes);
    std::vector<double> cell_ages(num_nodes);
    std::vector<double> cell_cycle_phases(num_nodes);
    std::vector<std::vector<double> > cellwise_data;

    if (CellwiseData<DIM>::Instance()->IsSetUp())
    {
        CellwiseData<DIM>* p_data = CellwiseData<DIM>::Instance();
        unsigned num_variables = p_data->GetNumVariables();
        for (unsigned var=0; var<num_variables; var++)
        {
            std::vector<double> cellwise_data_var(num_nodes);
            cellwise_data.push_back(cellwise_data_var);
        }
    }

    // Loop over cells
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        // Get the node index corresponding to this cell
        unsigned node_index = this->GetLocationIndexUsingCell(*cell_iter);

        if (this->mOutputCellAncestors)
        {
            double ancestor_index = (cell_iter->GetAncestor() == UNSIGNED_UNSET) ? (-1.0) : (double)cell_iter->GetAncestor();
            cell_ancestors[node_index] = ancestor_index;
        }
        if (this->mOutputCellProliferativeTypes)
        {
            double cell_type = cell_iter->GetCellCycleModel()->GetCellProliferativeType();
            cell_types[node_index] = cell_type;
        }
        if (this->mOutputCellMutationStates)
        {
            double mutation_state = cell_iter->GetMutationState()->GetColour();
            cell_mutation_states[node_index] = mutation_state;
        }
        if (this->mOutputCellAges)
        {
            double age = cell_iter->GetAge();
            cell_ages[node_index] = age;
        }
        if (this->mOutputCellCyclePhases)
        {
            double cycle_phase = cell_iter->GetCellCycleModel()->GetCurrentCellCyclePhase();
            cell_cycle_phases[node_index] = cycle_phase;
        }
        if (CellwiseData<DIM>::Instance()->IsSetUp())
        {
            CellwiseData<DIM>* p_data = CellwiseData<DIM>::Instance();
            unsigned num_variables = p_data->GetNumVariables();
            for (unsigned var=0; var<num_variables; var++)
            {
                cellwise_data[var][node_index] = p_data->GetValue(*cell_iter, var);
            }
        }
    }

    if (this->mOutputCellProliferativeTypes)
    {
        mesh_writer.AddPointData("Cell types", cell_types);
    }
    if (this->mOutputCellAncestors)
    {
        mesh_writer.AddPointData("Ancestors", cell_ancestors);
    }
    if (this->mOutputCellMutationStates)
    {
        mesh_writer.AddPointData("Mutation states", cell_mutation_states);
    }
    if (this->mOutputCellAges)
    {
        mesh_writer.AddPointData("Ages", cell_ages);
    }
    if (this->mOutputCellCyclePhases)
    {
        mesh_writer.AddPointData("Cycle phases", cell_cycle_phases);
    }
    if (this->mOutputCellVolumes)
    {
//        mesh_writer.AddCellData("Cell volumes", cell_volumes);
        ///\todo (#1598) ??
    }
    if (CellwiseData<DIM>::Instance()->IsSetUp())
    {
        for (unsigned var=0; var<cellwise_data.size(); var++)
        {
            std::stringstream data_name;
            data_name << "Cellwise data " << var;
            std::vector<double> cellwise_data_var = cellwise_data[var];
            mesh_writer.AddPointData(data_name.str(), cellwise_data_var);
        }
    }
    
    {
        //Make a copy of the nodes in a disposable mesh for writing...
        TetrahedralMesh<DIM,DIM> mesh;
        mesh.ConstructNodesWithoutMesh(mNodes);
        mesh_writer.WriteFilesUsingMesh(mesh);
    }
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

template class NodeBasedCellPopulation<1>;
template class NodeBasedCellPopulation<2>;
template class NodeBasedCellPopulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NodeBasedCellPopulation)
