#ifndef TISSUE_CPP
#define TISSUE_CPP

#include "Tissue.hpp"
#include "CancerParameters.hpp"
#include "VoronoiTessellation.cpp"

enum cell_colours
{
    STEM_COLOUR, // 0
    TRANSIT_COLOUR, // 1
    DIFFERENTIATED_COLOUR, // 2
    EARLY_CANCER_COLOUR, // 3
    LATE_CANCER_COLOUR, // 4
    LABELLED_COLOUR, // 5
    APOPTOSIS_COLOUR, // 6
    INVISIBLE_COLOUR, // visualizer treats '7' as invisible
    SPECIAL_LABEL_START
};

///\todo: Make this constructor take in ghost nodes, and validate the three objects
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

    // This must always be true    
    assert( mCells.size() <= mrMesh.GetNumNodes() );

    // Set up the node map
    for (std::list<TissueCell>::iterator it = mCells.begin();
         it != mCells.end();
         ++it)
    {
        /// \todo Check it points to a real cell; if not do
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


// Check every node either has a cell associated with it or is a ghost node
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
    // Find the node to which this cell corresponds
    unsigned node_index = rCell.GetNodeIndex();
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
    // Reinitialise all to false..
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
            // Check if this cell is in a marked spring
            std::vector<const std::set<TissueCell*>*> pairs_to_remove; // Pairs that must be purged
            for (std::set<std::set<TissueCell*> >::iterator it1 = mMarkedSprings.begin();
                 it1 != mMarkedSprings.end();
                 ++it1)
            {
                const std::set<TissueCell*>& r_pair = *it1;
                for (std::set<TissueCell*>::iterator it2 = r_pair.begin();
                     it2 != r_pair.end();
                     ++it2)
                {
                    TissueCell* p_cell = *it2;
                    if (p_cell == &(*it))
                    {
                        // Remember to purge this spring
                        pairs_to_remove.push_back(&r_pair);
                        break;
                    }
                }
            }
            // Purge any marked springs that contained this cell
            for (std::vector<const std::set<TissueCell*>* >::iterator pair_it = pairs_to_remove.begin();
                 pair_it != pairs_to_remove.end();
                 ++pair_it)
            {
                mMarkedSprings.erase(**pair_it);
            }
            
            // Remove the node from the mesh
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
    
    // There is reason not to substract one position from the other (cyclidrical meshes). clever gary
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
        // Copy mIsGhostNode. nodes bool
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
            
            // This shouldn't ever happen, as the cell vectors is only ever living cells
            assert(!map.IsDeleted(old_node_index));
            unsigned new_node_index = map.GetNewIndex(old_node_index);
            it->SetNodeIndex(new_node_index);
            mNodeCellMap[new_node_index] = &(*it);
        }
    }
    
    // Purge any marked springs that are no longer springs
    std::vector<const std::set<TissueCell*>*> springs_to_remove;
    for (std::set<std::set<TissueCell*> >::iterator spring_it = mMarkedSprings.begin();
         spring_it != mMarkedSprings.end();
         ++spring_it)
    {
        const std::set<TissueCell*>& r_pair = *spring_it;
        assert(r_pair.size() == 2);
        TissueCell* p_cell_1 = *(r_pair.begin());
        TissueCell* p_cell_2 = *(++r_pair.begin());
        Node<DIM>* p_node_1 = GetNodeCorrespondingToCell(*p_cell_1);
        Node<DIM>* p_node_2 = GetNodeCorrespondingToCell(*p_cell_2);
        
        bool joined = false;
        // For each element containing node1, if it also contains node2 then the cells are joined
        std::set<unsigned> node2_elements = p_node_2->rGetContainingElementIndices();
        for (typename Node<DIM>::ContainingElementIterator elt_it = p_node_1->ContainingElementsBegin();
             elt_it != p_node_1->ContainingElementsEnd();
             ++elt_it)
        {
            unsigned elt_index = *elt_it;
            if (node2_elements.find(elt_index) != node2_elements.end())
            {
                joined = true;
                break;
            }
        }
        
        // If no longer joined, remove this spring from the set
        if (!joined)
        {
            springs_to_remove.push_back(&r_pair);
        }
    }
    // Remove any springs necessary
    for (std::vector<const std::set<TissueCell*>* >::iterator spring_it = springs_to_remove.begin();
         spring_it != springs_to_remove.end();
         ++spring_it)
    {
        mMarkedSprings.erase(**spring_it);
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

template<unsigned DIM> 
void Tissue<DIM>::SetCellAncestorsToNodeIndices()
{
    for(Iterator cell_iter = Begin(); cell_iter!=End(); ++cell_iter)
    {
        cell_iter->SetAncestor(cell_iter->GetNodeIndex());
    }
}

template<unsigned DIM> 
void Tissue<DIM>::SetBottomCellAncestors()
{
    unsigned index = 0;
    for(Iterator cell_iter = Begin(); cell_iter!=End(); ++cell_iter)
    {
        if (cell_iter.rGetLocation()[1] < 0.5)
        {
            cell_iter->SetAncestor(index++);
        }
    }
}

template<unsigned DIM> 
std::set<unsigned> Tissue<DIM>::GetCellAncestors()
{
    std::set<unsigned> remaining_ancestors;
    for(Iterator cell_iter = Begin(); cell_iter!=End(); ++cell_iter)
    {
        remaining_ancestors.insert(cell_iter->GetAncestor());
    }
    return remaining_ancestors;
}


//////////////////////////////////////////////////////////////////////////////
//                             Iterator class                               // 
//////////////////////////////////////////////////////////////////////////////
template<unsigned DIM>
TissueCell& Tissue<DIM>::Iterator::operator*()
{
    assert(!IsAtEnd());
    return *mCellIter;
}

template<unsigned DIM>
TissueCell* Tissue<DIM>::Iterator::operator->()
{
    assert(!IsAtEnd());
    return &(*mCellIter);
}

template<unsigned DIM>
Node<DIM>* Tissue<DIM>::Iterator::GetNode()
{
    assert(!IsAtEnd());
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
        if (!IsAtEnd())
        {
            mNodeIndex = mCellIter->GetNodeIndex();
        }
    }
    while (!IsAtEnd() && !IsRealCell());
  
    return (*this);
}

template<unsigned DIM>
bool Tissue<DIM>::Iterator::IsRealCell()
{
    assert(mrTissue.rGetGhostNodes().size() == mrTissue.rGetMesh().GetNumAllNodes() );
    return !(mrTissue.rGetGhostNodes()[mNodeIndex] || GetNode()->IsDeleted() || (*this)->IsDead());
}

template<unsigned DIM>
bool Tissue<DIM>::Iterator::IsAtEnd()
{
    return mCellIter == mrTissue.rGetCells().end();
}

template<unsigned DIM>
Tissue<DIM>::Iterator::Iterator(Tissue& rTissue, std::list<TissueCell>::iterator cellIter)
    : mrTissue(rTissue),
      mCellIter(cellIter)
{
    // Make sure the tissue isn't empty
    assert(mrTissue.rGetCells().size() > 0);
    if (!IsAtEnd())
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
//                             Output methods                               // 
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM> 
c_vector<unsigned,5> Tissue<DIM>::GetCellTypeCount()
{
    return mCellTypeCount;
}

template<unsigned DIM>
void Tissue<DIM>::CreateOutputFiles(const std::string &rDirectory, bool rCleanOutputDirectory, bool outputCellTypes)
{
    OutputFileHandler output_file_handler(rDirectory, rCleanOutputDirectory);
    mpNodeFile = output_file_handler.OpenOutputFile("results.viznodes");
    mpElementFile = output_file_handler.OpenOutputFile("results.vizelements");
    mpCellTypesFile = output_file_handler.OpenOutputFile("celltypes.dat");
    
    if (outputCellTypes)
    {
        *mpCellTypesFile <<   "Time\t Healthy\t Labelled\t APC_1\t APC_2\t BETA_CAT \n";
    }
}

template<unsigned DIM>
void Tissue<DIM>::CloseOutputFiles()
{
    mpNodeFile->close();
    mpElementFile->close();
    mpCellTypesFile->close();
}

template<unsigned DIM>  
void Tissue<DIM>::WriteResultsToFiles(bool outputCellTypes)
{
    // Write current simulation time
    SimulationTime *p_simulation_time = SimulationTime::Instance();
    double time = p_simulation_time->GetDimensionalisedTime();
    unsigned cell_counter[5];
    for(unsigned i=0; i < 5; i++)
    {
        cell_counter[i] =0;
    }
    
    *mpNodeFile <<  time << "\t";
    *mpElementFile <<  time << "\t";
    
    if (outputCellTypes)
    {
        *mpCellTypesFile <<  time << "\t";
    }

    // Write node files
    for (unsigned index = 0; index<mrMesh.GetNumAllNodes(); index++)
    {
        unsigned colour = STEM_COLOUR; // all green if no cells have been passed in
         
        if (mIsGhostNode[index]==true)
        {
            colour = INVISIBLE_COLOUR; // visualizer treats '7' as invisible
        }
        else if (mrMesh.GetNode(index)->IsDeleted())
        {
            // Do nothing
        }
        else if (mNodeCellMap[index]->GetAncestor()!=UNSIGNED_UNSET)
        {
            TissueCell* p_cell = mNodeCellMap[index];
            colour = SPECIAL_LABEL_START + p_cell->GetAncestor();
        }
/// \todo remove this if - facade eventually shouldn't be able to have empty cells vector
        else if (mCells.size()>0)
        {
            TissueCell* p_cell = mNodeCellMap[index];
            
            CellType type = p_cell->GetCellType();
            CellMutationState mutation = p_cell->GetMutationState();
            
            // Set colours dependent on Stem, Transit, Differentiate, HepaOne
            if (type == STEM)
            {
                colour = STEM_COLOUR;
            }
            else if (type == TRANSIT)
            {
                colour = TRANSIT_COLOUR;
            }
            else
            {
                colour = DIFFERENTIATED_COLOUR;
            }
            
            // Override colours for mutant or labelled cells.
            if (mutation != HEALTHY && mutation != ALARCON_NORMAL)
            {
                if (mutation == LABELLED || mutation == ALARCON_CANCER)
                {
                    colour = LABELLED_COLOUR;
                    if (outputCellTypes)
                    {
                        cell_counter[1]++;
                    }
                }
                if (mutation == APC_ONE_HIT)
                {
                    colour = EARLY_CANCER_COLOUR;
                    if (outputCellTypes)
                    {
                        cell_counter[2]++;
                    }
                }
                if (mutation == APC_TWO_HIT )
                {
                    colour = LATE_CANCER_COLOUR;
                    if (outputCellTypes)
                    {
                        cell_counter[3]++;
                    }  
                }
                if ( mutation == BETA_CATENIN_ONE_HIT)
                {
                    colour = LATE_CANCER_COLOUR;
                    if (outputCellTypes)
                    {
                        cell_counter[4]++;
                    }  
                }
            }
            else // It's healthy, or normal in the sense of the Alarcon model
            {
                if (outputCellTypes)
                {
                    cell_counter[0]++;
                }  
            }
            
            if (p_cell->HasApoptosisBegun())
            {   
                // For any type of cell set the colour to this if it is undergoing apoptosis.
                colour = APOPTOSIS_COLOUR;   
            }
        }
        
        if (!mrMesh.GetNode(index)->IsDeleted())
        {
            const c_vector<double,DIM>& position = mrMesh.GetNode(index)->rGetLocation();
            
            for(unsigned i=0; i<DIM; i++)
            {
                *mpNodeFile << position[i] << " ";
            }
            *mpNodeFile << colour << " ";
        }
    }
    
    // Write element data files
    for (unsigned elem_index = 0; elem_index<mrMesh.GetNumAllElements(); elem_index++)
    {
        if (!mrMesh.GetElement(elem_index)->IsDeleted())
        {
            for(unsigned i=0; i<DIM+1; i++)
            {
                *mpElementFile << mrMesh.GetElement(elem_index)->GetNodeGlobalIndex(i)<< " ";
            }
        }
    }
        
    if (outputCellTypes)
    {
        for(unsigned i=0; i < 5; i++)
        {
            mCellTypeCount[i] = cell_counter[i];
            *mpCellTypesFile <<  cell_counter[i] << "\t";
        }
        *mpCellTypesFile <<  "\n";
    }
    
    *mpNodeFile << "\n";
    *mpElementFile << "\n";
}


//////////////////////////////////////////////////////////////////////////////
//                          Spring iterator class                           // 
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
        
        // Check cell exists in tissue
        unsigned node_index = p_cell->GetNodeIndex();
        std::cout << "Cell at node " << node_index << " addr " << p_cell << std::endl << std::flush;
        TissueCell& r_cell = rGetCellAtNodeIndex(node_index);
        if (&r_cell != p_cell)
        {
            std::cout << "  Mismatch with tissue" << std::endl << std::flush;
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
    
    res = true;
    for (std::set<std::set<TissueCell*> >::iterator it1 = mMarkedSprings.begin();
         it1 != mMarkedSprings.end();
         ++it1)
    {
        const std::set<TissueCell*>& r_pair = *it1;
        assert(r_pair.size() == 2);
        for (std::set<TissueCell*>::iterator it2 = r_pair.begin();
             it2 != r_pair.end();
             ++it2)
        {
            TissueCell* p_cell = *it2;
            assert(p_cell);
            AbstractCellCycleModel *p_model = p_cell->GetCellCycleModel();
            assert(p_model);
            unsigned node_index = p_cell->GetNodeIndex();
            std::cout << "Cell at node " << node_index << " addr " << p_cell << std::endl << std::flush;
            
            // Check cell is alive
            if (p_cell->IsDead())
            {
                std::cout << "  Cell is dead" << std::endl << std::flush;
                res = false;
            }
            
            // Check cell exists in tissue
            TissueCell& r_cell = rGetCellAtNodeIndex(node_index);
            if (&r_cell != p_cell)
            {
                std::cout << "  Mismatch with tissue" << std::endl << std::flush;
                res = false;
            }
            
            // Check model links back to cell
            if (p_model->GetCell() != p_cell)
            {
                std::cout << "  Mismatch with cycle model" << std::endl << std::flush;
                res = false;
            }
        }
    }
    assert(res);
}
#undef COVERAGE_IGNORE

template<unsigned DIM>
std::set<TissueCell*> Tissue<DIM>::CreateCellPair(TissueCell& rCell1, TissueCell& rCell2)
{
    std::set<TissueCell *> cell_pair;
    cell_pair.insert(&rCell1);
    cell_pair.insert(&rCell2);
    return cell_pair;
}

template<unsigned DIM>
bool Tissue<DIM>::IsMarkedSpring(TissueCell& rCell1, TissueCell& rCell2)
{
    std::set<TissueCell *> cell_pair = CreateCellPair(rCell1, rCell2);
    return mMarkedSprings.find(cell_pair) != mMarkedSprings.end();
}

template<unsigned DIM>
void Tissue<DIM>::MarkSpring(TissueCell& rCell1, TissueCell& rCell2)
{
    std::set<TissueCell *> cell_pair = CreateCellPair(rCell1, rCell2);
    mMarkedSprings.insert(cell_pair);
}

template<unsigned DIM>
void Tissue<DIM>::UnmarkSpring(TissueCell& rCell1, TissueCell& rCell2)
{
    std::set<TissueCell *> cell_pair = CreateCellPair(rCell1, rCell2);
    mMarkedSprings.erase(cell_pair);
}

#endif //TISSUE_CPP

