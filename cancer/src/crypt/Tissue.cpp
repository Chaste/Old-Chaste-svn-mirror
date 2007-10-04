#ifndef CRYPT_CPP
#define CRYPT_CPP

#include "Tissue.hpp"
#include "CancerParameters.hpp"
#include "VoronoiTessellation.cpp"

///\todo: make this constructor take in ghost nodes, and validate the three objects
// are in sync ie num cells + num ghost nodes = num_nodes ? this would mean all ghosts
// *cannot* be cells, making it more difficult to construct the cells.
// also check cell.GetNodeIndices() is in the mesh, and covers the mesh, etc.
template<unsigned DIM>
Tissue<DIM>::Tissue(ConformingTetrahedralMesh<DIM, DIM>& rMesh,
                  const std::vector<TissueCell>& rCells,
                  bool deleteMesh)
             : mrMesh(rMesh),
               mCells(rCells.begin(), rCells.end())
{
    mpVoronoiTessellation = NULL;
    
    mDeleteMesh = deleteMesh;
    mIsGhostNode = std::vector<bool>(mrMesh.GetNumNodes(), false);
    
    mMaxCells = 10*mrMesh.GetNumNodes();
    mMaxElements = 10*mrMesh.GetNumElements();

    // this must always be true    
    assert( mCells.size() <= mrMesh.GetNumNodes() );

    // Set up the node map
    for (std::list<TissueCell>::iterator it = mCells.begin();
         it != mCells.end();
         ++it)
    {
        /// \todo check it points to a real cell; if not do
        /// it = mCells.erase(it); --it; continue;
        unsigned node_index = it->GetNodeIndex();
        mNodeCellMap[node_index] = &(*it);
    }
    
    for(unsigned i=0; i < 5; i++)
    {
        mCellTypeCount[i] =0;
    }
        
	Validate();
}

template<unsigned DIM>
Tissue<DIM>::Tissue(ConformingTetrahedralMesh<DIM, DIM>& rMesh)
             : mrMesh(rMesh)
{
    mpVoronoiTessellation = NULL;
    mDeleteMesh = true;
}

template<unsigned DIM>
Tissue<DIM>::~Tissue()
{
    if (mDeleteMesh)
    {
        delete &mrMesh;
    }
    
    
    delete mpVoronoiTessellation;
}

template<unsigned DIM>
void Tissue<DIM>::InitialiseCells()
{
    for(std::list<TissueCell>::iterator iter = mCells.begin();
        iter != mCells.end();
        ++iter)
    {
        iter->InitialiseCellCycleModel();
    }
}


// check every node either has a cell associated with it or is a ghost node
// (for the time being, we are allowing ghost nodes to also have cells 
// associated with it, although this isn't very clean)
template<unsigned DIM>
void Tissue<DIM>::Validate()
{
	std::vector<bool> validated_node = mIsGhostNode; 
	for(Iterator cell_iter = Begin(); cell_iter!=End(); ++cell_iter)
	{
		unsigned node_index = cell_iter->GetNodeIndex();
		validated_node[node_index] = true;
	}
	
	for(unsigned i=0; i<validated_node.size(); i++)
	{
		if(!validated_node[i])
		{
			std::stringstream ss;
			ss << "Node " << i << " does not appear to be a ghost node or have a cell associated with it";
            EXCEPTION(ss.str()); 
		}
	}
}

template<unsigned DIM>
TissueCell& Tissue<DIM>::rGetCellAtNodeIndex(unsigned nodeGlobalIndex)
{
    return *(mNodeCellMap[nodeGlobalIndex]);
}

template<unsigned DIM>
Node<DIM>* Tissue<DIM>::GetNodeCorrespondingToCell(const TissueCell& rCell)
{
    // find the node to which this cell corresponds
    std::map<unsigned, TissueCell*>::iterator it=mNodeCellMap.begin();
    while (it != mNodeCellMap.end() && (*it).second != &rCell)
    {
        it++;
    }
    assert (it != mNodeCellMap.end());
    unsigned node_index = (*it).first;
    return mrMesh.GetNode(node_index);   
}


template<unsigned DIM>
c_vector<double, DIM> Tissue<DIM>::GetLocationOfCell(const TissueCell& rCell)
{
    return GetNodeCorrespondingToCell(rCell)->rGetLocation();
}



template<unsigned DIM>
ConformingTetrahedralMesh<DIM, DIM>& Tissue<DIM>::rGetMesh()
{
    return mrMesh;
}

template<unsigned DIM>
std::list<TissueCell>& Tissue<DIM>::rGetCells()
{
    return mCells;
}

template<unsigned DIM>
const ConformingTetrahedralMesh<DIM, DIM>& Tissue<DIM>::rGetMesh() const
{
    return mrMesh;
}

template<unsigned DIM>
const std::list<TissueCell>& Tissue<DIM>::rGetCells() const
{
    return mCells;
}

template<unsigned DIM>
std::vector<bool>& Tissue<DIM>::rGetGhostNodes()
{
    return mIsGhostNode;
}

template<unsigned DIM>
std::set<unsigned> Tissue<DIM>::GetGhostNodeIndices()
{
    std::set<unsigned> ghost_node_indices;
    for (unsigned i=0; i<mIsGhostNode.size(); i++)
    {
        if (mIsGhostNode[i])
        {
            ghost_node_indices.insert(i);    
        }        
    }
    return ghost_node_indices;        
}

template<unsigned DIM>
void Tissue<DIM>::SetGhostNodes(const std::vector<bool>& rGhostNodes)
{
    mIsGhostNode = rGhostNodes;
}

template<unsigned DIM> 
void Tissue<DIM>::SetGhostNodes(const std::set<unsigned>& ghostNodeIndices)
{
    // reinitialise all to false..
    mIsGhostNode = std::vector<bool>(mrMesh.GetNumNodes(), false);

    // ..then update which ones are ghosts
    std::set<unsigned>::iterator iter = ghostNodeIndices.begin();
    while(iter!=ghostNodeIndices.end())
    {
        mIsGhostNode[*iter]=true;
        iter++;
    }
}


template<unsigned DIM>
unsigned Tissue<DIM>::RemoveDeadCells()
{
    unsigned num_removed = 0;
    for (std::list<TissueCell>::iterator it = mCells.begin();
         it != mCells.end();
         ++it)
    {
        if (it->IsDead())
        {
            num_removed++;
            mrMesh.DeleteNodePriorToReMesh(it->GetNodeIndex());
            it = mCells.erase(it);
            --it;
        }
    }
    return num_removed;
}

template<unsigned DIM>
void Tissue<DIM>::UpdateGhostPositions(double dt)
{
    std::vector<c_vector<double, DIM> > drdt(mrMesh.GetNumAllNodes());
    for (unsigned i=0; i<drdt.size(); i++)
    {
        drdt[i]=zero_vector<double>(DIM);
    }

    for(typename ConformingTetrahedralMesh<DIM, DIM>::EdgeIterator edge_iterator=mrMesh.EdgesBegin();
        edge_iterator!=mrMesh.EdgesEnd();
        ++edge_iterator)
    {
        unsigned nodeA_global_index = edge_iterator.GetNodeA()->GetIndex();
        unsigned nodeB_global_index = edge_iterator.GetNodeB()->GetIndex();

        c_vector<double, DIM> force = CalculateForceBetweenNodes(nodeA_global_index,nodeB_global_index);
         
        double damping_constant = CancerParameters::Instance()->GetDampingConstantNormal();
        
        
        if (!mIsGhostNode[nodeA_global_index])
        {
            drdt[nodeB_global_index] -= force / damping_constant;
        }
        else
        {
            drdt[nodeA_global_index] += force / damping_constant;
                
            if (mIsGhostNode[nodeB_global_index])
            {
                drdt[nodeB_global_index] -= force / damping_constant;
            }
        }
    }
    
    for (unsigned index = 0; index<mrMesh.GetNumAllNodes(); index++)
    {
        if ((!mrMesh.GetNode(index)->IsDeleted()) && mIsGhostNode[index])
        {
            ChastePoint<DIM> new_point(mrMesh.GetNode(index)->rGetLocation() + dt*drdt[index]);
            mrMesh.SetNode(index, new_point, false);
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
c_vector<double, DIM> Tissue<DIM>::CalculateForceBetweenNodes(const unsigned& rNodeAGlobalIndex, const unsigned& rNodeBGlobalIndex)
{
    assert(rNodeAGlobalIndex!=rNodeBGlobalIndex);
    c_vector<double, DIM> unit_difference;
    c_vector<double, DIM> node_a_location = mrMesh.GetNode(rNodeAGlobalIndex)->rGetLocation();
    c_vector<double, DIM> node_b_location = mrMesh.GetNode(rNodeBGlobalIndex)->rGetLocation();
    
    // there is reason not to substract one position from the other (cyclidrical meshes). clever gary
    unit_difference = mrMesh.GetVectorFromAtoB(node_a_location, node_b_location);   
    
    double distance_between_nodes = norm_2(unit_difference);
    
    unit_difference /= distance_between_nodes;
    
    double rest_length = 1.0;
    
    return CancerParameters::Instance()->GetSpringStiffness() * unit_difference * (distance_between_nodes - rest_length);
}

template<unsigned DIM>
void Tissue<DIM>::MoveCell(Iterator iter, ChastePoint<DIM>& rNewLocation)
{
    unsigned index = iter.GetNode()->GetIndex();
    mrMesh.SetNode(index, rNewLocation, false);
}

template<unsigned DIM>  
TissueCell* Tissue<DIM>::AddCell(TissueCell newCell, c_vector<double,DIM> newLocation)
{
    Node<DIM>* p_new_node = new Node<DIM>(mrMesh.GetNumNodes(), newLocation, false);   // never on boundary
              
    unsigned new_node_index = mrMesh.AddNode(p_new_node);

    newCell.SetNodeIndex(new_node_index);
    mCells.push_back(newCell);
    
    TissueCell *p_created_cell=&(mCells.back());
    mNodeCellMap[new_node_index] = p_created_cell;
    
    // Update size of IsGhostNode if necessary
    if (mrMesh.GetNumNodes() > mIsGhostNode.size())
    {
        mIsGhostNode.resize(mrMesh.GetNumNodes());
        mIsGhostNode[new_node_index] = false;
    }
    return p_created_cell;
}


template<unsigned DIM>
void Tissue<DIM>::ReMesh()
{
    NodeMap map(mrMesh.GetNumAllNodes());
    mrMesh.ReMesh(map);

    if(!map.IsIdentityMap())
    {
        // copy mIsGhostNode. nodes bool
        std::vector<bool> ghost_nodes_before_remesh = mIsGhostNode;    
        mIsGhostNode.clear();
        mIsGhostNode.resize(mrMesh.GetNumNodes());
        
        for(unsigned old_index=0; old_index<map.Size(); old_index++)
        {
            if(!map.IsDeleted(old_index))
            {
                unsigned new_index = map.GetNewIndex(old_index);
                mIsGhostNode[new_index] = ghost_nodes_before_remesh[old_index];
            }
        }

        // Fix up the mappings between cells and nodes
        mNodeCellMap.clear();
        for (std::list<TissueCell>::iterator it = mCells.begin();
             it != mCells.end();
             ++it)
        {
            unsigned old_node_index = it->GetNodeIndex();
            // this shouldn't ever happen, as the cell vectors is only ever living cells
            assert(!map.IsDeleted(old_node_index));
            unsigned new_node_index = map.GetNewIndex(old_node_index);
            it->SetNodeIndex(new_node_index);
            mNodeCellMap[new_node_index] = &(*it);
        }
    }
    
    Validate();
}

template<unsigned DIM>
unsigned Tissue<DIM>::GetNumRealCells()
{
	unsigned counter = 0;
	for(Iterator cell_iter = Begin(); cell_iter!=End(); ++cell_iter)
	{
		counter++;
	}
	return counter;
}


//////////////////////////////////////////////////////////////////////////////
//                             iterator class                               // 
//////////////////////////////////////////////////////////////////////////////
template<unsigned DIM>
TissueCell& Tissue<DIM>::Iterator::operator*()
{
    assert((*this) != mrTissue.End());
    return *mCellIter;
}


template<unsigned DIM>
TissueCell* Tissue<DIM>::Iterator::operator->()
{
    assert((*this) != mrTissue.End());
    return &(*mCellIter);
}


template<unsigned DIM>
Node<DIM>* Tissue<DIM>::Iterator::GetNode()
{
    assert((*this) != mrTissue.End());
    return mrTissue.rGetMesh().GetNode(mNodeIndex);
}

template<unsigned DIM>
const c_vector<double, DIM>& Tissue<DIM>::Iterator::rGetLocation()
{
    return GetNode()->rGetLocation();
}

template<unsigned DIM>
bool Tissue<DIM>::Iterator::operator!=(const Tissue<DIM>::Iterator& other)
{
    return mCellIter != other.mCellIter;   
}

template<unsigned DIM>
typename Tissue<DIM>::Iterator& Tissue<DIM>::Iterator::operator++()
{
    do
    {
        ++mCellIter;
        if((*this) != mrTissue.End())
        {
            mNodeIndex = mCellIter->GetNodeIndex();
        }
    }
    while ((*this) != mrTissue.End() && !IsRealCell());
  
    return (*this);
}

template<unsigned DIM>
bool Tissue<DIM>::Iterator::IsRealCell()
{
    assert(mrTissue.rGetGhostNodes().size() == mrTissue.rGetMesh().GetNumAllNodes() );
    return !(mrTissue.rGetGhostNodes()[mNodeIndex] || GetNode()->IsDeleted() || (*this)->IsDead());
}

template<unsigned DIM>
Tissue<DIM>::Iterator::Iterator(Tissue& rTissue, std::list<TissueCell>::iterator cellIter)
    : mrTissue(rTissue),
      mCellIter(cellIter)
{
    // Make sure the crypt isn't empty
    assert(mrTissue.rGetCells().size() > 0);
    if (mCellIter != mrTissue.rGetCells().end())
    {
        mNodeIndex = cellIter->GetNodeIndex();
    }
    // Make sure we start at a real cell
    if (mCellIter == mrTissue.rGetCells().begin() && !IsRealCell())
    {
        ++(*this);
    }
}

template<unsigned DIM>
typename Tissue<DIM>::Iterator Tissue<DIM>::Begin()
{
    return Iterator(*this, mCells.begin());
}

template<unsigned DIM>
typename Tissue<DIM>::Iterator Tissue<DIM>::End()
{
    return Iterator(*this, mCells.end());
}


//////////////////////////////////////////////////////////////////////////////
//                             output methods                               // 
//////////////////////////////////////////////////////////////////////////////
template<unsigned DIM>  
void Tissue<DIM>::SetMaxCells(unsigned maxCells)
{
    mMaxCells = maxCells;
    if (maxCells<mrMesh.GetNumAllNodes())
    {
        #define COVERAGE_IGNORE
        EXCEPTION("mMaxCells is less than the number of cells in the mesh.");
        #undef COVERAGE_IGNORE
    }
}

/**
 * Sets the maximum number of elements that the simulation will contain (for use by the datawriter)
 * default value is set to 10x the initial mesh value by the constructor.
 */
template<unsigned DIM> 
void Tissue<DIM>::SetMaxElements(unsigned maxElements)
{
    mMaxElements = maxElements;
    if (maxElements<mrMesh.GetNumAllElements())
    {
        #define COVERAGE_IGNORE
        EXCEPTION("mMaxElements is less than the number of elements in the mesh.");
        #undef COVERAGE_IGNORE
    }
}

template<unsigned DIM>  
void Tissue<DIM>::SetupTabulatedWriters(ColumnDataWriter& rNodeWriter, ColumnDataWriter& rElementWriter)
{   
    // set up node writer
    mNodeVarIds.time = rNodeWriter.DefineUnlimitedDimension("Time","hours");
    
    mNodeVarIds.types.resize(mMaxCells);
    mNodeVarIds.position_id.resize(mMaxCells);
    
    // set up per-cell variables
    for (unsigned cell=0; cell<mMaxCells; cell++)
    {
        std::stringstream cell_type_var_name, cell_x_position_var_name, cell_y_position_var_name, cell_z_position_var_name;
        cell_type_var_name << "cell_type_" << cell;
        
        cell_x_position_var_name << "cell_x_position_" << cell;
        cell_y_position_var_name << "cell_y_position_" << cell;
        cell_z_position_var_name << "cell_z_position_" << cell; // not used in 2d

        std::vector<std::string> cell_position_var_name_string;
        cell_position_var_name_string.push_back(cell_x_position_var_name.str());
        cell_position_var_name_string.push_back(cell_y_position_var_name.str());
        cell_position_var_name_string.push_back(cell_z_position_var_name.str()); // not used in 2d
        
        mNodeVarIds.types[cell]=rNodeWriter.DefineVariable(cell_type_var_name.str(),"dimensionless");
        for(unsigned i=0; i<DIM; i++)
        {
            mNodeVarIds.position_id[cell](i)=rNodeWriter.DefineVariable(cell_position_var_name_string[i],"rest_spring_length");
        }
    }
    
    rNodeWriter.EndDefineMode();

    // set up element writer
    mElemVarIds.time = rElementWriter.DefineUnlimitedDimension("Time","hours");
    
    // Set up columns for element writer
    mElemVarIds.node_id.resize(mMaxElements);
    
    for (unsigned elem_index = 0; elem_index<mMaxElements; elem_index++)
    {
        std::stringstream nodeA_var_name, nodeB_var_name, nodeC_var_name, nodeD_var_name;
        
        nodeA_var_name << "nodeA_" << elem_index;
        nodeB_var_name << "nodeB_" << elem_index;
        nodeC_var_name << "nodeC_" << elem_index;
        nodeD_var_name << "nodeD_" << elem_index;

        std::vector<std::string> node_var_name_string;
        node_var_name_string.push_back(nodeA_var_name.str());
        node_var_name_string.push_back(nodeB_var_name.str());
        node_var_name_string.push_back(nodeC_var_name.str());
        node_var_name_string.push_back(nodeD_var_name.str());
        
        for(unsigned i=0; i<DIM+1; i++)
        {
            mElemVarIds.node_id[elem_index](i) = rElementWriter.DefineVariable(node_var_name_string[i],"dimensionless");
        }
    }
    
    rElementWriter.EndDefineMode();
}


template<unsigned DIM> 
c_vector<unsigned,5> Tissue<DIM>::GetCellTypeCount()
{
    return mCellTypeCount;
}

template<unsigned DIM>  
void Tissue<DIM>::WriteResultsToFiles(ColumnDataWriter& rNodeWriter,
                                     ColumnDataWriter& rElementWriter,
                                     std::ofstream& rNodeFile, 
                                     std::ofstream& rElementFile,
                                     std::ofstream& rCellTypesFile,
                                     bool writeTabulatedResults,
                                     bool writeVisualizerResults,
                                     bool OutputCellTypes)
{
    // Write current simulation time
    SimulationTime *p_simulation_time = SimulationTime::Instance();
    double time = p_simulation_time->GetDimensionalisedTime();
    unsigned cell_counter[5];
    for(unsigned i=0; i < 5; i++)
    {
        cell_counter[i] =0;
    }
    
    if (writeVisualizerResults)
    {
        rNodeFile <<  time << "\t";
        rElementFile <<  time << "\t";
    }
    
    if (OutputCellTypes)
    {
        rCellTypesFile <<  time << "\t";
    }
       
    if (writeTabulatedResults)
    {
        rNodeWriter.PutVariable(mNodeVarIds.time, time);
        rElementWriter.PutVariable(mElemVarIds.time, time);
    }
    
        
    // write node files
    for (unsigned index = 0; index<mrMesh.GetNumAllNodes(); index++)
    {
        if (index>mMaxCells)
        {
            #define COVERAGE_IGNORE
            EXCEPTION("\nNumber of cells exceeds mMaxCells. Use SetMaxCells(unsigned) to increase it.\n");
            #undef COVERAGE_IGNORE
        }
        unsigned colour = 0; // all green if no cells have been passed in
        
        if (mIsGhostNode[index]==true)
        {
            colour = 7; // visualizer treats '7' as invisible
        }
        else if (mrMesh.GetNode(index)->IsDeleted())
        {
            // do nothing
        }
/// \todo remove this if - facade eventually shouldn't be able to have empty cells vector
        else if (mCells.size()>0)
        {
            TissueCell* p_cell = mNodeCellMap[index];
            
            CryptCellType type = p_cell->GetCellType();
            CryptCellMutationState mutation = p_cell->GetMutationState();
            
            // Set colours dependent on Stem, Transit, Differentiate
            if (type == STEM)
            {
                colour = 0;
            }
            else if (type == TRANSIT)
            {
                colour = 1;
            }
            else
            {
                colour = 2;
            }
            
            // Override colours for mutant or labelled cells.
            if (mutation != HEALTHY)
            {
                if (mutation == LABELLED)
                {
                    colour = 5;
                    if (OutputCellTypes)
                    {
                        cell_counter[1]++;
                    }
                }
                if (mutation == APC_ONE_HIT)
                {
                    colour = 3;
                    if (OutputCellTypes)
                    {
                        cell_counter[2]++;
                    }
                }
                if (mutation == APC_TWO_HIT )
                {
                    colour = 4;
                    if (OutputCellTypes)
                    {
                        cell_counter[3]++;
                    }  
                }
                if ( mutation == BETA_CATENIN_ONE_HIT)
                {
                    colour = 4;
                    if (OutputCellTypes)
                    {
                        cell_counter[4]++;
                    }  
                }
            }
            else // its healthy
            {
                if (OutputCellTypes)
                {
                    cell_counter[0]++;
                }  
            }
            
            if (p_cell->HasApoptosisBegun())
            {   // For any type of cell set the colour to this if it is undergoing apoptosis.
                colour = 6;   
            }
            
            
        }
        
        if (!mrMesh.GetNode(index)->IsDeleted())
        {
            const c_vector<double,DIM>& position = mrMesh.GetNode(index)->rGetLocation();
            if (writeVisualizerResults)
            {
                for(unsigned i=0; i<DIM; i++)
                {
                    rNodeFile << position[i] << " ";
                }
                rNodeFile << colour << " ";
            }
            if (writeTabulatedResults)
            {
                for(unsigned i=0; i<DIM; i++)
                {
                    rNodeWriter.PutVariable(mNodeVarIds.position_id[index](i), position[i]);
                }
                rNodeWriter.PutVariable(mNodeVarIds.types[index], colour);
            }
        }
    }
    
    // write element data files
    for (unsigned elem_index = 0; elem_index<mrMesh.GetNumAllElements(); elem_index++)
    {
        if (elem_index>mMaxElements)
        {
            #define COVERAGE_IGNORE
            EXCEPTION("Maximum number of elements (mMaxElements) exceeded.\nUse SetMaxElements(unsigned) to increase it.\n");
            #undef COVERAGE_IGNORE
        }
        if (!mrMesh.GetElement(elem_index)->IsDeleted())
        {
            if (writeVisualizerResults)
            {
                for(unsigned i=0; i<DIM+1; i++)
                {
                    rElementFile << mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(i)<< " ";
                }
            }
            if (writeTabulatedResults)
            {
                for(unsigned i=0; i<DIM+1; i++)
                {
                    rElementWriter.PutVariable(mElemVarIds.node_id[elem_index](i), mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(i));
                }
            }
        }
    }
    
    
    if (OutputCellTypes)
    {
        for(unsigned i=0; i < 5; i++)
        {
            mCellTypeCount[i] = cell_counter[i];
            rCellTypesFile <<  cell_counter[i] << "\t";
        }
        rCellTypesFile <<  "\n";
    }
    
    if (writeVisualizerResults)
    {
        rNodeFile << "\n";
        rElementFile << "\n";
    }
    if (writeTabulatedResults)
    {
        rNodeWriter.AdvanceAlongUnlimitedDimension();
        rElementWriter.AdvanceAlongUnlimitedDimension();
    }
}


//////////////////////////////////////////////////////////////////////////////
//                          spring iterator class                           // 
//////////////////////////////////////////////////////////////////////////////
template<unsigned DIM>
Node<DIM>* Tissue<DIM>::SpringIterator::GetNodeA()
{
    return mEdgeIter.GetNodeA();
}

template<unsigned DIM>
Node<DIM>* Tissue<DIM>::SpringIterator::GetNodeB()
{
    return mEdgeIter.GetNodeB();
}

template<unsigned DIM>
TissueCell& Tissue<DIM>::SpringIterator::rGetCellA()
{
    assert((*this) != mrTissue.SpringsEnd());
    return mrTissue.rGetCellAtNodeIndex(mEdgeIter.GetNodeA()->GetIndex());
}


template<unsigned DIM>
TissueCell& Tissue<DIM>::SpringIterator::rGetCellB()
{
    assert((*this) != mrTissue.SpringsEnd());
    return mrTissue.rGetCellAtNodeIndex(mEdgeIter.GetNodeB()->GetIndex());
}


template<unsigned DIM>
bool Tissue<DIM>::SpringIterator::operator!=(const Tissue<DIM>::SpringIterator& other)
{
    return (mEdgeIter != other.mEdgeIter);
}

template<unsigned DIM>
typename Tissue<DIM>::SpringIterator& Tissue<DIM>::SpringIterator::operator++()
{
    bool edge_is_ghost = false;
    
    do
    {
        ++mEdgeIter;
        if(*this !=mrTissue.SpringsEnd())
        {
            bool a_is_ghost = mrTissue.mIsGhostNode[mEdgeIter.GetNodeA()->GetIndex()];
            bool b_is_ghost = mrTissue.mIsGhostNode[mEdgeIter.GetNodeB()->GetIndex()];

            edge_is_ghost = (a_is_ghost || b_is_ghost);
        }
    }
    while( *this!=mrTissue.SpringsEnd() && edge_is_ghost ); 

    return (*this);
}


template<unsigned DIM>
Tissue<DIM>::SpringIterator::SpringIterator(Tissue& rTissue,
                                           typename ConformingTetrahedralMesh<DIM,DIM>::EdgeIterator edgeIter)
    : mrTissue(rTissue),
      mEdgeIter(edgeIter)
{
    if(mEdgeIter!=mrTissue.mrMesh.EdgesEnd())
    {
        bool a_is_ghost = mrTissue.mIsGhostNode[mEdgeIter.GetNodeA()->GetIndex()];
        bool b_is_ghost = mrTissue.mIsGhostNode[mEdgeIter.GetNodeB()->GetIndex()];

        if(a_is_ghost || b_is_ghost)
        {
            ++(*this);
        }
    }
}

template<unsigned DIM>
typename Tissue<DIM>::SpringIterator Tissue<DIM>::SpringsBegin()
{
    return SpringIterator(*this, mrMesh.EdgesBegin());
}

template<unsigned DIM>
typename Tissue<DIM>::SpringIterator Tissue<DIM>::SpringsEnd()
{
    return SpringIterator(*this, mrMesh.EdgesEnd());
}

template<unsigned DIM>
void Tissue<DIM>::CreateVoronoiTessellation()
{
    delete mpVoronoiTessellation;
    mpVoronoiTessellation = new VoronoiTessellation<DIM>(mrMesh);
}

template<unsigned DIM>
VoronoiTessellation<DIM>& Tissue<DIM>::rGetVoronoiTessellation()
{
    assert(mpVoronoiTessellation!=NULL);
    return *mpVoronoiTessellation;
}

#define COVERAGE_IGNORE
template<unsigned DIM>
void Tissue<DIM>::CheckTissueCellPointers()
{
    bool res=true;
    for (std::list<TissueCell>::iterator it=mCells.begin();
        it!=mCells.end();
        ++it)
    {
        TissueCell* p_cell=&(*it);
        assert(p_cell);
        AbstractCellCycleModel *p_model = p_cell->GetCellCycleModel();
        assert(p_model);
        // Check cell exists in crypt
        unsigned node_index = p_cell->GetNodeIndex();
        std::cout << "Cell at node " << node_index << " addr " << p_cell << std::endl << std::flush;
        TissueCell& r_cell = rGetCellAtNodeIndex(node_index);
        if (&r_cell != p_cell)
        {
            std::cout << "  Mismatch with crypt" << std::endl << std::flush;
            res = false;
        }
        // Check model links back to cell
        if (p_model->GetCell() != p_cell)
        {
            std::cout << "  Mismatch with cycle model" << std::endl << std::flush;
            res = false;
        }
    }
    assert(res);
}
#undef COVERAGE_IGNORE

#endif //CRYPT_CPP

