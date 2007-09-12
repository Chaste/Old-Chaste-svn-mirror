#ifndef CRYPT_CPP
#define CRYPT_CPP

#include "Crypt.hpp"
#include "CancerParameters.hpp"

///\todo: make this constructor take in ghost nodes, and validate the three objects
// are in sync ie num cells + num ghost nodes = num_nodes ? this would mean all ghosts
// *cannot* be cells, making it more difficult to construct the cells.
// also check cell.GetNodeIndices() is in the mesh, and covers the mesh, etc.
template<unsigned DIM>
Crypt<DIM>::Crypt(ConformingTetrahedralMesh<DIM, DIM>& rMesh,
                  const std::vector<MeinekeCryptCell>& rCells,
                  bool deleteMesh)
             : mrMesh(rMesh),
               mCells(rCells.begin(), rCells.end())
{
    mDeleteMesh = deleteMesh;
    mIsGhostNode = std::vector<bool>(mrMesh.GetNumNodes(), false);
    
    mMaxCells = 10*mrMesh.GetNumNodes();
    mMaxElements = 10*mrMesh.GetNumElements();

    // this must always be true    
    assert( mCells.size() <= mrMesh.GetNumNodes() );

    // Set up the node map
    for (std::list<MeinekeCryptCell>::iterator it = mCells.begin();
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
        
    // remove this line when facade is finished 
    //   - no don't, validate does important stuff like check each non-ghost 
    //     node has a cell
	Validate();
}

template<unsigned DIM>
Crypt<DIM>::~Crypt()
{
    if (mDeleteMesh)
    {
        delete &mrMesh;
    }
}


// check every node either has a cell associated with it or is a ghost node
// (for the time being, we are allowing ghost nodes to also have cells 
// associated with it, although this isn't very clean)
template<unsigned DIM>
void Crypt<DIM>::Validate()
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
MeinekeCryptCell& Crypt<DIM>::rGetCellAtNodeIndex(unsigned nodeGlobalIndex)
{
    return *(mNodeCellMap[nodeGlobalIndex]);
}

template<unsigned DIM>
c_vector<double, DIM> Crypt<DIM>::GetLocationOfCell(const MeinekeCryptCell& rCell)
{
    // find the node to which this cell corresponds
    std::map<unsigned, MeinekeCryptCell*>::iterator it=mNodeCellMap.begin();
    while (it != mNodeCellMap.end() && (*it).second != &rCell)
    {
        it++;
    }
    assert (it != mNodeCellMap.end());
    unsigned node_index = (*it).first;
    return mrMesh.GetNode(node_index)->rGetLocation();   
}

template<unsigned DIM>
ConformingTetrahedralMesh<DIM, DIM>& Crypt<DIM>::rGetMesh()
{
    return mrMesh;
}

template<unsigned DIM>
std::list<MeinekeCryptCell>& Crypt<DIM>::rGetCells()
{
    return mCells;
}

template<unsigned DIM>
const ConformingTetrahedralMesh<DIM, DIM>& Crypt<DIM>::rGetMesh() const
{
    return mrMesh;
}

template<unsigned DIM>
const std::list<MeinekeCryptCell>& Crypt<DIM>::rGetCells() const
{
    return mCells;
}

template<unsigned DIM>
std::vector<bool>& Crypt<DIM>::rGetGhostNodes()
{
    return mIsGhostNode;
}

template<unsigned DIM>
std::set<unsigned> Crypt<DIM>::GetGhostNodeIndices()
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
void Crypt<DIM>::SetGhostNodes(const std::vector<bool>& rGhostNodes)
{
    mIsGhostNode = rGhostNodes;
}

template<unsigned DIM> 
void Crypt<DIM>::SetGhostNodes(const std::set<unsigned>& ghostNodeIndices)
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
unsigned Crypt<DIM>::RemoveDeadCells()
{
    unsigned num_removed = 0;
    for (std::list<MeinekeCryptCell>::iterator it = mCells.begin();
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
void Crypt<DIM>::UpdateGhostPositions(double dt)
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
c_vector<double, DIM> Crypt<DIM>::CalculateForceBetweenNodes(const unsigned& rNodeAGlobalIndex, const unsigned& rNodeBGlobalIndex)
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
void Crypt<DIM>::MoveCell(Iterator iter, ChastePoint<DIM>& rNewLocation)
{
    unsigned index = iter.GetNode()->GetIndex();
    mrMesh.SetNode(index, rNewLocation, false);
}

template<unsigned DIM>  
MeinekeCryptCell* Crypt<DIM>::AddCell(MeinekeCryptCell newCell, c_vector<double,DIM> newLocation)
{
    Node<DIM>* p_new_node = new Node<DIM>(mrMesh.GetNumNodes(), newLocation, false);   // never on boundary
              
    unsigned new_node_index = mrMesh.AddNode(p_new_node);

    newCell.SetNodeIndex(new_node_index);
    mCells.push_back(newCell);
    
    MeinekeCryptCell *p_created_cell=&(mCells.back());
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
void Crypt<DIM>::ReMesh()
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
        for (std::list<MeinekeCryptCell>::iterator it = mCells.begin();
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
unsigned Crypt<DIM>::GetNumRealCells()
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
MeinekeCryptCell& Crypt<DIM>::Iterator::operator*()
{
    assert((*this) != mrCrypt.End());
    return *mCellIter;
}


template<unsigned DIM>
MeinekeCryptCell* Crypt<DIM>::Iterator::operator->()
{
    assert((*this) != mrCrypt.End());
    return &(*mCellIter);
}


template<unsigned DIM>
Node<DIM>* Crypt<DIM>::Iterator::GetNode()
{
    assert((*this) != mrCrypt.End());
    return mrCrypt.rGetMesh().GetNode(mNodeIndex);
}

template<unsigned DIM>
const c_vector<double, DIM>& Crypt<DIM>::Iterator::rGetLocation()
{
    return GetNode()->rGetLocation();
}

template<unsigned DIM>
bool Crypt<DIM>::Iterator::operator!=(const Crypt<DIM>::Iterator& other)
{
    return mCellIter != other.mCellIter;   
}

template<unsigned DIM>
typename Crypt<DIM>::Iterator& Crypt<DIM>::Iterator::operator++()
{
    do
    {
        ++mCellIter;
        if((*this) != mrCrypt.End())
        {
            mNodeIndex = mCellIter->GetNodeIndex();
        }
    }
    while ((*this) != mrCrypt.End() && !IsRealCell());
  
    return (*this);
}

template<unsigned DIM>
bool Crypt<DIM>::Iterator::IsRealCell()
{
    assert(mrCrypt.rGetGhostNodes().size() == mrCrypt.rGetMesh().GetNumAllNodes() );
    return !(mrCrypt.rGetGhostNodes()[mNodeIndex] || GetNode()->IsDeleted() || (*this)->IsDead());
}

template<unsigned DIM>
Crypt<DIM>::Iterator::Iterator(Crypt& rCrypt, std::list<MeinekeCryptCell>::iterator cellIter)
    : mrCrypt(rCrypt),
      mCellIter(cellIter)
{
    // Make sure the crypt isn't empty
    assert(mrCrypt.rGetCells().size() > 0);
    if (mCellIter != mrCrypt.rGetCells().end())
    {
        mNodeIndex = cellIter->GetNodeIndex();
    }
    // Make sure we start at a real cell
    if (mCellIter == mrCrypt.rGetCells().begin() && !IsRealCell())
    {
        ++(*this);
    }
}

template<unsigned DIM>
typename Crypt<DIM>::Iterator Crypt<DIM>::Begin()
{
    return Iterator(*this, mCells.begin());
}

template<unsigned DIM>
typename Crypt<DIM>::Iterator Crypt<DIM>::End()
{
    return Iterator(*this, mCells.end());
}


//////////////////////////////////////////////////////////////////////////////
//                             output methods                               // 
//////////////////////////////////////////////////////////////////////////////
template<unsigned DIM>  
void Crypt<DIM>::SetMaxCells(unsigned maxCells)
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
void Crypt<DIM>::SetMaxElements(unsigned maxElements)
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
void Crypt<DIM>::SetupTabulatedWriters(ColumnDataWriter& rNodeWriter, ColumnDataWriter& rElementWriter)
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
c_vector<unsigned,5> Crypt<DIM>::GetCellTypeCount()
{
    return mCellTypeCount;
}

template<unsigned DIM>  
void Crypt<DIM>::WriteResultsToFiles(ColumnDataWriter& rNodeWriter,
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
            MeinekeCryptCell* p_cell = mNodeCellMap[index];
            
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
Node<DIM>* Crypt<DIM>::SpringIterator::GetNodeA()
{
    return mEdgeIter.GetNodeA();
}

template<unsigned DIM>
Node<DIM>* Crypt<DIM>::SpringIterator::GetNodeB()
{
    return mEdgeIter.GetNodeB();
}

template<unsigned DIM>
MeinekeCryptCell& Crypt<DIM>::SpringIterator::rGetCellA()
{
    assert((*this) != mrCrypt.SpringsEnd());
    return mrCrypt.rGetCellAtNodeIndex(mEdgeIter.GetNodeA()->GetIndex());
}


template<unsigned DIM>
MeinekeCryptCell& Crypt<DIM>::SpringIterator::rGetCellB()
{
    assert((*this) != mrCrypt.SpringsEnd());
    return mrCrypt.rGetCellAtNodeIndex(mEdgeIter.GetNodeB()->GetIndex());
}


template<unsigned DIM>
bool Crypt<DIM>::SpringIterator::operator!=(const Crypt<DIM>::SpringIterator& other)
{
    return (mEdgeIter != other.mEdgeIter);
}

template<unsigned DIM>
typename Crypt<DIM>::SpringIterator& Crypt<DIM>::SpringIterator::operator++()
{
    bool edge_is_ghost = false;
    
    do
    {
        ++mEdgeIter;
        if(*this !=mrCrypt.SpringsEnd())
        {
            bool a_is_ghost = mrCrypt.mIsGhostNode[mEdgeIter.GetNodeA()->GetIndex()];
            bool b_is_ghost = mrCrypt.mIsGhostNode[mEdgeIter.GetNodeB()->GetIndex()];

            edge_is_ghost = (a_is_ghost || b_is_ghost);
        }
    }
    while( *this!=mrCrypt.SpringsEnd() && edge_is_ghost ); 

    return (*this);
}


template<unsigned DIM>
Crypt<DIM>::SpringIterator::SpringIterator(Crypt& rCrypt,
                                           typename ConformingTetrahedralMesh<DIM,DIM>::EdgeIterator edgeIter)
    : mrCrypt(rCrypt),
      mEdgeIter(edgeIter)
{
    if(mEdgeIter!=mrCrypt.mrMesh.EdgesEnd())
    {
        bool a_is_ghost = mrCrypt.mIsGhostNode[mEdgeIter.GetNodeA()->GetIndex()];
        bool b_is_ghost = mrCrypt.mIsGhostNode[mEdgeIter.GetNodeB()->GetIndex()];

        if(a_is_ghost || b_is_ghost)
        {
            ++(*this);
        }
    }
}

template<unsigned DIM>
typename Crypt<DIM>::SpringIterator Crypt<DIM>::SpringsBegin()
{
    return SpringIterator(*this, mrMesh.EdgesBegin());
}

template<unsigned DIM>
typename Crypt<DIM>::SpringIterator Crypt<DIM>::SpringsEnd()
{
    return SpringIterator(*this, mrMesh.EdgesEnd());
}


#endif //CRYPT_CPP

