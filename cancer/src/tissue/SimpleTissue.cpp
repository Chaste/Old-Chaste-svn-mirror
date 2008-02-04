#ifndef SIMPLETISSUE_CPP
#define SIMPLETISSUE_CPP

#include "SimpleTissue.hpp"

template<unsigned DIM>
SimpleTissue<DIM>::SimpleTissue(const std::vector<Node<DIM> >& rNodes, 
                                const std::vector<TissueCell>& rCells)
        : AbstractTissue<DIM>(rCells),
          mNodes(rNodes)
{
    // Set up the node map
    for (std::list<TissueCell>::iterator it = this->mCells.begin();
         it != this->mCells.end();
         ++it)
    {
        // \todo Check it points to a real cell; if not do
        // it = this->mCells.erase(it); --it; continue;
        unsigned node_index = it->GetNodeIndex();
        this->mNodeCellMap[node_index] = &(*it);
    }
    
    for(unsigned i=0; i<5; i++)
    {
        this->mCellTypeCount[i]=0;
    }
    
    Validate();    
}

template<unsigned DIM>
void SimpleTissue<DIM>::Validate()
{
    std::vector<bool> validated_node(GetNumNodes()); 
    
    for (Iterator cell_iter=Begin(); cell_iter!=End(); ++cell_iter)
    {
        unsigned node_index = cell_iter->GetNodeIndex();
        validated_node[node_index] = true;
    }
    
    for (unsigned i=0; i<validated_node.size(); i++)
    {
        if(!validated_node[i])
        {
            std::stringstream ss;
            ss << "Node " << i << " does not appear to have a cell associated with it";
            EXCEPTION(ss.str()); 
        }
    }
}

template<unsigned DIM>
void SimpleTissue<DIM>::InitialiseCells()
{
    for(std::list<TissueCell>::iterator iter = this->mCells.begin();
        iter != this->mCells.end();
        ++iter)
    {
        iter->InitialiseCellCycleModel();
    }
}

template<unsigned DIM>
std::vector<Node<DIM> >& SimpleTissue<DIM>::rGetNodes()
{
    return mNodes;
}

template<unsigned DIM>
Node<DIM>* SimpleTissue<DIM>::GetNode(unsigned index)
{
    return &(mNodes[index]);    
}

template<unsigned DIM>
Node<DIM>* SimpleTissue<DIM>::GetNodeCorrespondingToCell(const TissueCell& rCell)
{
    // Find the node to which this cell corresponds
    unsigned node_index = rCell.GetNodeIndex();
    return GetNode(node_index);   
}

template<unsigned DIM>
c_vector<double, DIM> SimpleTissue<DIM>::GetLocationOfCell(const TissueCell& rCell)
{
    return GetNodeCorrespondingToCell(rCell)->rGetLocation();
}

template<unsigned DIM>
void SimpleTissue<DIM>::SetNode(unsigned index, ChastePoint<DIM> point)
{
    mNodes[index].SetPoint(point);    
}

template<unsigned DIM>
void SimpleTissue<DIM>::MoveCell(Iterator iter, ChastePoint<DIM>& rNewLocation)
{
    unsigned index = iter.GetNode()->GetIndex();
    SetNode(index, rNewLocation);
}

template<unsigned DIM>
TissueCell* SimpleTissue<DIM>::AddCell(TissueCell newCell, c_vector<double,DIM> newLocation)
{
    // Create a new node
    Node<DIM>* p_new_node = new Node<DIM>(GetNumNodes(), newLocation, false); // never on boundary
    
    unsigned new_node_index = AddNode(p_new_node);

    // Associate the new cell with the node
    newCell.SetNodeIndex(new_node_index);
    this->mCells.push_back(newCell);
    
    TissueCell *p_created_cell = &(this->mCells.back());
    this->mNodeCellMap[new_node_index] = p_created_cell;
        
    return p_created_cell;
}

template<unsigned DIM>
void SimpleTissue<DIM>::RemoveNode(unsigned index)
{
    mNodes.erase(mNodes.begin() + index);    
}

template<unsigned DIM>
unsigned SimpleTissue<DIM>::RemoveDeadCells()
{
    unsigned num_removed = 0;
    
    for (std::list<TissueCell>::iterator it = this->mCells.begin();
         it != this->mCells.end();
         ++it)
    {
        if (it->IsDead())
        {   
            // Remove the node from the mesh
            num_removed++;
            RemoveNode(it->GetNodeIndex());
            it = this->mCells.erase(it);
            --it;
        }
    }
    return num_removed;
}

template<unsigned DIM>
TissueCell& SimpleTissue<DIM>::rGetCellAtNodeIndex(unsigned index)
{
    return *(this->mNodeCellMap[index]);
}

template<unsigned DIM>
unsigned SimpleTissue<DIM>::AddNode(Node<DIM> *pNewNode)
{
    // \todo: employ a std::vector of deleted node indices to re-use indices? 
    pNewNode->SetIndex(mNodes.size());
    mNodes.push_back(*pNewNode);        
    return pNewNode->GetIndex();
}    

template<unsigned DIM>
void SimpleTissue<DIM>::UpdateNodeCellMap()
{
    // We expect the node indices to be {0,1,...,num_nodes}
    std::set<unsigned> expected_node_indices;    
    for (unsigned i=0; i<GetNumNodes(); i++)
    {
        expected_node_indices.insert(i);
    }
    
    // Get the actual set of node indices
    std::set<unsigned> node_indices;    
    for (Iterator cell_iter=Begin(); cell_iter!=End(); ++cell_iter)
    {
        unsigned node_index = cell_iter->GetNodeIndex();
        node_indices.insert(node_index);       
    }
    
    // If necessary, update the node cell map        
    if (node_indices != expected_node_indices)
    {
        // Fix up the mappings between cells and nodes
        this->mNodeCellMap.clear();
        unsigned new_node_index = 0;
        for (Iterator cell_iter=Begin(); cell_iter!=End(); ++cell_iter)
        {
            cell_iter->SetNodeIndex(new_node_index);
            this->mNodeCellMap[new_node_index] = &(*cell_iter);
            new_node_index++;
        }
    }
}

template<unsigned DIM>
unsigned SimpleTissue<DIM>::GetNumRealCells()
{
    unsigned counter = 0;
    for(Iterator cell_iter=Begin(); cell_iter!=End(); ++cell_iter)
    {
        counter++;
    }
    return counter;
}

template<unsigned DIM>
unsigned SimpleTissue<DIM>::GetNumNodes()
{
    return mNodes.size();
}

template<unsigned DIM>
void SimpleTissue<DIM>::SetCellAncestorsToNodeIndices()
{
    for(Iterator cell_iter = Begin(); cell_iter!=End(); ++cell_iter)
    {
        cell_iter->SetAncestor(cell_iter->GetNodeIndex());
    }
}

template<unsigned DIM>
std::set<unsigned> SimpleTissue<DIM>::GetCellAncestors()
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
TissueCell& SimpleTissue<DIM>::Iterator::operator*()
{
    assert(!IsAtEnd());
    return *mCellIter;
}

template<unsigned DIM>
TissueCell* SimpleTissue<DIM>::Iterator::operator->()
{
    assert(!IsAtEnd());
    return &(*mCellIter);
}

template<unsigned DIM>
Node<DIM>* SimpleTissue<DIM>::Iterator::GetNode()
{
    assert(!IsAtEnd());
    return mrTissue.GetNode(mNodeIndex);
}

template<unsigned DIM>
const c_vector<double, DIM>& SimpleTissue<DIM>::Iterator::rGetLocation()
{
    return GetNode()->rGetLocation();
}

template<unsigned DIM>
bool SimpleTissue<DIM>::Iterator::operator!=(const SimpleTissue<DIM>::Iterator& other)
{
    return mCellIter != other.mCellIter;   
}

template<unsigned DIM>
typename SimpleTissue<DIM>::Iterator& SimpleTissue<DIM>::Iterator::operator++()
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
bool SimpleTissue<DIM>::Iterator::IsRealCell()
{
    return !(GetNode()->IsDeleted() || (*this)->IsDead());
}

template<unsigned DIM>
bool SimpleTissue<DIM>::Iterator::IsAtEnd()
{
    return mCellIter == mrTissue.rGetCells().end();
}

template<unsigned DIM>
SimpleTissue<DIM>::Iterator::Iterator(SimpleTissue& rTissue, std::list<TissueCell>::iterator cellIter)
    : mrTissue(rTissue),
      mCellIter(cellIter)
{
    // Make sure the tissue isn't empty
    assert(mrTissue.rGetCells().size() > 0);
    if (!IsAtEnd())
    {
        mNodeIndex = cellIter->GetNodeIndex();
    }
}


template<unsigned DIM>
typename SimpleTissue<DIM>::Iterator SimpleTissue<DIM>::Begin()
{
    return Iterator(*this, this->mCells.begin());
}

template<unsigned DIM>
typename SimpleTissue<DIM>::Iterator SimpleTissue<DIM>::End()
{
    return Iterator(*this, this->mCells.end());
}


//////////////////////////////////////////////////////////////////////////////
//                             Output methods                               // 
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void SimpleTissue<DIM>::CreateOutputFiles(const std::string &rDirectory, bool rCleanOutputDirectory, bool outputCellTypes)
{
    OutputFileHandler output_file_handler(rDirectory, rCleanOutputDirectory);
    this->mpNodeFile = output_file_handler.OpenOutputFile("results.viznodes");
    this->mpCellTypesFile = output_file_handler.OpenOutputFile("celltypes.dat");
    
    if (outputCellTypes)
    {
        *this->mpCellTypesFile <<   "Time\t Healthy\t Labelled\t APC_1\t APC_2\t BETA_CAT \n";
    }
}

template<unsigned DIM>
void SimpleTissue<DIM>::CloseOutputFiles()
{
    this->mpNodeFile->close();
    this->mpCellTypesFile->close();
}

template<unsigned DIM>  
void SimpleTissue<DIM>::WriteResultsToFiles(bool outputCellTypes)
{
    // Write current simulation time
    SimulationTime *p_simulation_time = SimulationTime::Instance();
    double time = p_simulation_time->GetDimensionalisedTime();
    unsigned cell_counter[5];
    for(unsigned i=0; i < 5; i++)
    {
        cell_counter[i] = 0;
    }
    
    *this->mpNodeFile <<  time << "\t";
    
    if (outputCellTypes)
    {
        *this->mpCellTypesFile <<  time << "\t";
    }

    // Write node files
    for (unsigned index=0; index<GetNumNodes(); index++)
    {
        unsigned colour = STEM_COLOUR; // all green if no cells have been passed in
         
        if (GetNode(index)->IsDeleted())
        {
            // Do nothing
        }
        else if (this->mNodeCellMap[index]->GetAncestor()!=UNSIGNED_UNSET)
        {
            TissueCell* p_cell = this->mNodeCellMap[index];
            colour = SPECIAL_LABEL_START + p_cell->GetAncestor();
        }
        else if (this->mCells.size() > 0)
        {
            TissueCell* p_cell = this->mNodeCellMap[index];
            
            CellType type = p_cell->GetCellType();
            CellMutationState mutation = p_cell->GetMutationState();
            
            // Set colours dependent on Stem, Transit, Differentiated, HepaOne
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
        if ( !(GetNode(index)->IsDeleted()) )
        {
            const c_vector<double,DIM>& position = GetNode(index)->rGetLocation();
            
            for(unsigned i=0; i<DIM; i++)
            {
                *this->mpNodeFile << position[i] << " ";
            }
            *this->mpNodeFile << colour << " ";
        }
    }
   
    if (outputCellTypes)
    {
        for(unsigned i=0; i < 5; i++)
        {
            this->mCellTypeCount[i] = cell_counter[i];
            *this->mpCellTypesFile <<  cell_counter[i] << "\t";
        }
        *this->mpCellTypesFile <<  "\n";
    }
    
    *this->mpNodeFile << "\n";
}


#endif //SIMPLETISSUE_CPP

#ifndef SIMPLETISSUE_CPP
#define SIMPLETISSUE_CPP

#include "SimpleTissue.hpp"

template<unsigned DIM>
SimpleTissue<DIM>::SimpleTissue(const std::vector<Node<DIM> >& rNodes, 
                                const std::vector<TissueCell>& rCells)
        : AbstractTissue<DIM>(rCells),
          mNodes(rNodes)
{
    // Set up the node map
    for (std::list<TissueCell>::iterator it = this->mCells.begin();
         it != this->mCells.end();
         ++it)
    {
        // \todo Check it points to a real cell; if not do
        // it = this->mCells.erase(it); --it; continue;
        unsigned node_index = it->GetNodeIndex();
        this->mNodeCellMap[node_index] = &(*it);
    }
    
    for(unsigned i=0; i<5; i++)
    {
        this->mCellTypeCount[i]=0;
    }
    
    Validate();    
}

template<unsigned DIM>
void SimpleTissue<DIM>::Validate()
{
    std::vector<bool> validated_node(GetNumNodes()); 
    
    for (Iterator cell_iter=Begin(); cell_iter!=End(); ++cell_iter)
    {
        unsigned node_index = cell_iter->GetNodeIndex();
        validated_node[node_index] = true;
    }
    
    for (unsigned i=0; i<validated_node.size(); i++)
    {
        if(!validated_node[i])
        {
            std::stringstream ss;
            ss << "Node " << i << " does not appear to have a cell associated with it";
            EXCEPTION(ss.str()); 
        }
    }
}

template<unsigned DIM>
void SimpleTissue<DIM>::InitialiseCells()
{
    for(std::list<TissueCell>::iterator iter = this->mCells.begin();
        iter != this->mCells.end();
        ++iter)
    {
        iter->InitialiseCellCycleModel();
    }
}

template<unsigned DIM>
std::vector<Node<DIM> >& SimpleTissue<DIM>::rGetNodes()
{
    return mNodes;
}

template<unsigned DIM>
Node<DIM>* SimpleTissue<DIM>::GetNode(unsigned index)
{
    return &(mNodes[index]);    
}

template<unsigned DIM>
Node<DIM>* SimpleTissue<DIM>::GetNodeCorrespondingToCell(const TissueCell& rCell)
{
    // Find the node to which this cell corresponds
    unsigned node_index = rCell.GetNodeIndex();
    return GetNode(node_index);   
}

template<unsigned DIM>
c_vector<double, DIM> SimpleTissue<DIM>::GetLocationOfCell(const TissueCell& rCell)
{
    return GetNodeCorrespondingToCell(rCell)->rGetLocation();
}

template<unsigned DIM>
void SimpleTissue<DIM>::SetNode(unsigned index, ChastePoint<DIM> point)
{
    mNodes[index].SetPoint(point);    
}

template<unsigned DIM>
void SimpleTissue<DIM>::MoveCell(Iterator iter, ChastePoint<DIM>& rNewLocation)
{
    unsigned index = iter.GetNode()->GetIndex();
    SetNode(index, rNewLocation);
}

template<unsigned DIM>
TissueCell* SimpleTissue<DIM>::AddCell(TissueCell newCell, c_vector<double,DIM> newLocation)
{
    // Create a new node
    Node<DIM>* p_new_node = new Node<DIM>(GetNumNodes(), newLocation, false); // never on boundary
    
    unsigned new_node_index = AddNode(p_new_node);

    // Associate the new cell with the node
    newCell.SetNodeIndex(new_node_index);
    this->mCells.push_back(newCell);
    
    TissueCell *p_created_cell = &(this->mCells.back());
    this->mNodeCellMap[new_node_index] = p_created_cell;
        
    return p_created_cell;
}

template<unsigned DIM>
void SimpleTissue<DIM>::RemoveNode(unsigned index)
{
    mNodes.erase(mNodes.begin() + index);    
}

template<unsigned DIM>
unsigned SimpleTissue<DIM>::RemoveDeadCells()
{
    unsigned num_removed = 0;
    
    for (std::list<TissueCell>::iterator it = this->mCells.begin();
         it != this->mCells.end();
         ++it)
    {
        if (it->IsDead())
        {   
            // Remove the node from the mesh
            num_removed++;
            RemoveNode(it->GetNodeIndex());
            it = this->mCells.erase(it);
            --it;
        }
    }
    return num_removed;
}

template<unsigned DIM>
TissueCell& SimpleTissue<DIM>::rGetCellAtNodeIndex(unsigned index)
{
    return *(this->mNodeCellMap[index]);
}

template<unsigned DIM>
unsigned SimpleTissue<DIM>::AddNode(Node<DIM> *pNewNode)
{
    // \todo: employ a std::vector of deleted node indices to re-use indices? 
    pNewNode->SetIndex(mNodes.size());
    mNodes.push_back(*pNewNode);        
    return pNewNode->GetIndex();
}    

template<unsigned DIM>
void SimpleTissue<DIM>::UpdateNodeCellMap()
{
    // We expect the node indices to be {0,1,...,num_nodes}
    std::set<unsigned> expected_node_indices;    
    for (unsigned i=0; i<GetNumNodes(); i++)
    {
        expected_node_indices.insert(i);
    }
    
    // Get the actual set of node indices
    std::set<unsigned> node_indices;    
    for (Iterator cell_iter=Begin(); cell_iter!=End(); ++cell_iter)
    {
        unsigned node_index = cell_iter->GetNodeIndex();
        node_indices.insert(node_index);       
    }
    
    // If necessary, update the node cell map        
    if (node_indices != expected_node_indices)
    {
        // Fix up the mappings between cells and nodes
        this->mNodeCellMap.clear();
        unsigned new_node_index = 0;
        for (Iterator cell_iter=Begin(); cell_iter!=End(); ++cell_iter)
        {
            cell_iter->SetNodeIndex(new_node_index);
            this->mNodeCellMap[new_node_index] = &(*cell_iter);
            new_node_index++;
        }
    }
}

template<unsigned DIM>
unsigned SimpleTissue<DIM>::GetNumRealCells()
{
    unsigned counter = 0;
    for(Iterator cell_iter=Begin(); cell_iter!=End(); ++cell_iter)
    {
        counter++;
    }
    return counter;
}

template<unsigned DIM>
unsigned SimpleTissue<DIM>::GetNumNodes()
{
    return mNodes.size();
}

template<unsigned DIM>
void SimpleTissue<DIM>::SetCellAncestorsToNodeIndices()
{
    for(Iterator cell_iter = Begin(); cell_iter!=End(); ++cell_iter)
    {
        cell_iter->SetAncestor(cell_iter->GetNodeIndex());
    }
}

template<unsigned DIM>
std::set<unsigned> SimpleTissue<DIM>::GetCellAncestors()
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
TissueCell& SimpleTissue<DIM>::Iterator::operator*()
{
    assert(!IsAtEnd());
    return *mCellIter;
}

template<unsigned DIM>
TissueCell* SimpleTissue<DIM>::Iterator::operator->()
{
    assert(!IsAtEnd());
    return &(*mCellIter);
}

template<unsigned DIM>
Node<DIM>* SimpleTissue<DIM>::Iterator::GetNode()
{
    assert(!IsAtEnd());
    return mrTissue.GetNode(mNodeIndex);
}

template<unsigned DIM>
const c_vector<double, DIM>& SimpleTissue<DIM>::Iterator::rGetLocation()
{
    return GetNode()->rGetLocation();
}

template<unsigned DIM>
bool SimpleTissue<DIM>::Iterator::operator!=(const SimpleTissue<DIM>::Iterator& other)
{
    return mCellIter != other.mCellIter;   
}

template<unsigned DIM>
typename SimpleTissue<DIM>::Iterator& SimpleTissue<DIM>::Iterator::operator++()
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
bool SimpleTissue<DIM>::Iterator::IsRealCell()
{
    return !(GetNode()->IsDeleted() || (*this)->IsDead());
}

template<unsigned DIM>
bool SimpleTissue<DIM>::Iterator::IsAtEnd()
{
    return mCellIter == mrTissue.rGetCells().end();
}

template<unsigned DIM>
SimpleTissue<DIM>::Iterator::Iterator(SimpleTissue& rTissue, std::list<TissueCell>::iterator cellIter)
    : mrTissue(rTissue),
      mCellIter(cellIter)
{
    // Make sure the tissue isn't empty
    assert(mrTissue.rGetCells().size() > 0);
    if (!IsAtEnd())
    {
        mNodeIndex = cellIter->GetNodeIndex();
    }
}


template<unsigned DIM>
typename SimpleTissue<DIM>::Iterator SimpleTissue<DIM>::Begin()
{
    return Iterator(*this, this->mCells.begin());
}

template<unsigned DIM>
typename SimpleTissue<DIM>::Iterator SimpleTissue<DIM>::End()
{
    return Iterator(*this, this->mCells.end());
}


//////////////////////////////////////////////////////////////////////////////
//                             Output methods                               // 
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void SimpleTissue<DIM>::CreateOutputFiles(const std::string &rDirectory, bool rCleanOutputDirectory, bool outputCellTypes)
{
    OutputFileHandler output_file_handler(rDirectory, rCleanOutputDirectory);
    this->mpNodeFile = output_file_handler.OpenOutputFile("results.viznodes");
    this->mpCellTypesFile = output_file_handler.OpenOutputFile("celltypes.dat");
    
    if (outputCellTypes)
    {
        *this->mpCellTypesFile <<   "Time\t Healthy\t Labelled\t APC_1\t APC_2\t BETA_CAT \n";
    }
}

template<unsigned DIM>
void SimpleTissue<DIM>::CloseOutputFiles()
{
    this->mpNodeFile->close();
    this->mpCellTypesFile->close();
}

template<unsigned DIM>  
void SimpleTissue<DIM>::WriteResultsToFiles(bool outputCellTypes)
{
    // Write current simulation time
    SimulationTime *p_simulation_time = SimulationTime::Instance();
    double time = p_simulation_time->GetDimensionalisedTime();
    unsigned cell_counter[5];
    for(unsigned i=0; i < 5; i++)
    {
        cell_counter[i] = 0;
    }
    
    *this->mpNodeFile <<  time << "\t";
    
    if (outputCellTypes)
    {
        *this->mpCellTypesFile <<  time << "\t";
    }

    // Write node files
    for (unsigned index=0; index<GetNumNodes(); index++)
    {
        unsigned colour = STEM_COLOUR; // all green if no cells have been passed in
         
        if (GetNode(index)->IsDeleted())
        {
            // Do nothing
        }
        else if (this->mNodeCellMap[index]->GetAncestor()!=UNSIGNED_UNSET)
        {
            TissueCell* p_cell = this->mNodeCellMap[index];
            colour = SPECIAL_LABEL_START + p_cell->GetAncestor();
        }
        else if (this->mCells.size() > 0)
        {
            TissueCell* p_cell = this->mNodeCellMap[index];
            
            CellType type = p_cell->GetCellType();
            CellMutationState mutation = p_cell->GetMutationState();
            
            // Set colours dependent on Stem, Transit, Differentiated, HepaOne
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
        if ( !(GetNode(index)->IsDeleted()) )
        {
            const c_vector<double,DIM>& position = GetNode(index)->rGetLocation();
            
            for(unsigned i=0; i<DIM; i++)
            {
                *this->mpNodeFile << position[i] << " ";
            }
            *this->mpNodeFile << colour << " ";
        }
    }
   
    if (outputCellTypes)
    {
        for(unsigned i=0; i < 5; i++)
        {
            this->mCellTypeCount[i] = cell_counter[i];
            *this->mpCellTypesFile <<  cell_counter[i] << "\t";
        }
        *this->mpCellTypesFile <<  "\n";
    }
    
    *this->mpNodeFile << "\n";
}


#endif //SIMPLETISSUE_CPP

