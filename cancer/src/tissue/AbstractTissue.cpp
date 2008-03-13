#ifndef ABSTRACTTISSUE_CPP
#define ABSTRACTTISSUE_CPP

#include "AbstractTissue.hpp"
#include "CancerParameters.hpp"

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

template<unsigned DIM>
AbstractTissue<DIM>::AbstractTissue(const std::vector<TissueCell>& rCells)
             : mCells(rCells.begin(), rCells.end()),               
               mTissueContainsMesh(false),
               mTissueContainsGhostNodes(false)
{
    // There must be at least one cell
    assert(mCells.size() > 0);
    
    // Set up the node map
    for (std::list<TissueCell>::iterator it = mCells.begin();
         it != mCells.end();
         ++it)
    {
        /// \todo Check it points to a real cell; if not do
        /// it = this->mCells.erase(it); --it; continue;
        unsigned node_index = it->GetNodeIndex();
        mNodeCellMap[node_index] = &(*it);
    }
    
    // Initialise cell counts to zero    
    for(unsigned i=0; i<NUM_CELL_MUTATION_STATES; i++)
    {
        mCellMutationStateCount[i] = 0;
    }
    for(unsigned i=0; i<NUM_CELL_TYPES; i++)
    {
        mCellTypeCount[i] = 0;
    }
    for(unsigned i=0; i<5; i++)
    {
        mCellCyclePhaseCount[i] = 0;
    }
}

template<unsigned DIM>
void AbstractTissue<DIM>::InitialiseCells()
{
    for (std::list<TissueCell>::iterator iter = mCells.begin();
        iter != mCells.end();
        ++iter)
    {
        iter->InitialiseCellCycleModel();
    }
}

template<unsigned DIM>
std::list<TissueCell>& AbstractTissue<DIM>::rGetCells()
{
    return mCells;
}

template<unsigned DIM>
const std::list<TissueCell>& AbstractTissue<DIM>::rGetCells() const
{
    return this->mCells;
}

template<unsigned DIM>
bool AbstractTissue<DIM>::HasMesh()
{
    return mTissueContainsMesh;
}

template<unsigned DIM>
bool AbstractTissue<DIM>::HasGhostNodes()
{
    return mTissueContainsGhostNodes;
}

template<unsigned DIM>
unsigned AbstractTissue<DIM>::GetNumRealCells()
{
    unsigned counter = 0;
    for (typename AbstractTissue<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        counter++;
    }
    return counter;
}

template<unsigned DIM>
Node<DIM>* AbstractTissue<DIM>::GetNodeCorrespondingToCell(const TissueCell& rCell)
{
    unsigned node_index = rCell.GetNodeIndex();
    return GetNode(node_index);
}

template<unsigned DIM>
c_vector<double, DIM> AbstractTissue<DIM>::GetLocationOfCell(const TissueCell& rCell)
{
    return GetNodeCorrespondingToCell(rCell)->rGetLocation();
}

template<unsigned DIM>
void AbstractTissue<DIM>::SetCellAncestorsToNodeIndices()
{
    for(typename AbstractTissue<DIM>::Iterator cell_iter = this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        cell_iter->SetAncestor(cell_iter->GetNodeIndex());
    }
}

template<unsigned DIM> 
std::set<unsigned> AbstractTissue<DIM>::GetCellAncestors()
{
    std::set<unsigned> remaining_ancestors;
    for (typename AbstractTissue<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        remaining_ancestors.insert(cell_iter->GetAncestor());
    }
    return remaining_ancestors;
}

template<unsigned DIM> 
c_vector<unsigned, NUM_CELL_MUTATION_STATES> AbstractTissue<DIM>::GetCellMutationStateCount()
{
    return mCellMutationStateCount;
}

template<unsigned DIM> 
c_vector<unsigned, NUM_CELL_TYPES> AbstractTissue<DIM>::GetCellTypeCount()
{
    return mCellTypeCount;
}

template<unsigned DIM> 
c_vector<unsigned, 5> AbstractTissue<DIM>::GetCellCyclePhaseCount()
{
    return mCellCyclePhaseCount;
}

template<unsigned DIM> 
unsigned AbstractTissue<DIM>::GetGhostNodesSize()
{
    return GetNumNodes();
}

template<unsigned DIM> 
bool AbstractTissue<DIM>::IsGhostNode(unsigned index)
{
    return false;
}    

template<unsigned DIM>
TissueCell& AbstractTissue<DIM>::rGetCellAtNodeIndex(unsigned index)
{
    return *(mNodeCellMap[index]);
}

//////////////////////////////////////////////////////////////////////////////
//                             Iterator class                               // 
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
TissueCell& AbstractTissue<DIM>::Iterator::operator*()
{
    assert(!IsAtEnd());
    return *mCellIter;
}

template<unsigned DIM>
TissueCell* AbstractTissue<DIM>::Iterator::operator->()
{
    assert(!IsAtEnd());
    return &(*mCellIter);
}

template<unsigned DIM>
Node<DIM>* AbstractTissue<DIM>::Iterator::GetNode()
{
    assert(!IsAtEnd());
    return mrTissue.GetNode(mNodeIndex);
}

template<unsigned DIM>
const c_vector<double, DIM>& AbstractTissue<DIM>::Iterator::rGetLocation()
{
    return GetNode()->rGetLocation();
}

template<unsigned DIM>
bool AbstractTissue<DIM>::Iterator::operator!=(const AbstractTissue<DIM>::Iterator& other)
{
    return mCellIter != other.mCellIter;   
}

template<unsigned DIM>
typename AbstractTissue<DIM>::Iterator& AbstractTissue<DIM>::Iterator::operator++()
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
bool AbstractTissue<DIM>::Iterator::IsRealCell()
{
    // TODO: move this assertion elsewhere, after mIsGhostNode has been deserialized! 
    // assert(mrTissue.GetGhostNodesSize() == mrTissue.GetNumNodes());
    return !(mrTissue.IsGhostNode(mNodeIndex) || GetNode()->IsDeleted() || (*this)->IsDead());
}

template<unsigned DIM>
bool AbstractTissue<DIM>::Iterator::IsAtEnd()
{
    return mCellIter == mrTissue.rGetCells().end();
}

template<unsigned DIM>
AbstractTissue<DIM>::Iterator::Iterator(AbstractTissue& rTissue, std::list<TissueCell>::iterator cellIter)
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
typename AbstractTissue<DIM>::Iterator AbstractTissue<DIM>::Begin()
{
    return Iterator(*this, this->mCells.begin());
}

template<unsigned DIM>
typename AbstractTissue<DIM>::Iterator AbstractTissue<DIM>::End()
{
    return Iterator(*this, this->mCells.end());
}

//////////////////////////////////////////////////////////////////////////////
//                             Output methods                               // 
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void AbstractTissue<DIM>::CreateOutputFiles(const std::string &rDirectory, 
                                            bool rCleanOutputDirectory, 
                                            bool outputCellMutationStates,
                                            bool outputCellTypes,
                                            bool outputCellVariables,
                                            bool outputCellCyclePhases)
{
    OutputFileHandler output_file_handler(rDirectory, rCleanOutputDirectory);
    mpNodeFile = output_file_handler.OpenOutputFile("results.viznodes");
    
    if (outputCellMutationStates)
    {
        mpCellMutationStatesFile = output_file_handler.OpenOutputFile("cellmutationstates.dat");
        *mpCellMutationStatesFile <<   "Time\t Healthy\t Labelled\t APC_1\t APC_2\t BETA_CAT \n";
    }
    if (outputCellTypes)
    {
        mpCellTypesFile = output_file_handler.OpenOutputFile("celltypes.dat");
    }
    if (outputCellVariables)
    {
        mpCellVariablesFile = output_file_handler.OpenOutputFile("cellvariables.dat");
    }
    if (outputCellCyclePhases)
    {
        mpCellCyclePhasesFile = output_file_handler.OpenOutputFile("cellcyclephases.dat");
    }
}

template<unsigned DIM>
void AbstractTissue<DIM>::CloseOutputFiles(bool outputCellMutationStates,
                                           bool outputCellTypes,
                                           bool outputCellVariables,
                                           bool outputCellCyclePhases)
{
    mpNodeFile->close();
    if (outputCellMutationStates)
    {
        mpCellMutationStatesFile->close();
    }
    if (outputCellTypes)
    {
        mpCellTypesFile->close();
    }
    if (outputCellVariables)
    {
        mpCellVariablesFile->close();
    }
    if (outputCellCyclePhases)
    {
        mpCellCyclePhasesFile->close();
    }
}

template<unsigned DIM>  
void AbstractTissue<DIM>::WriteResultsToFiles(bool outputCellMutationStates, 
                                              bool outputCellTypes, 
                                              bool outputCellVariables, 
                                              bool outputCellCyclePhases)
{   
    // Write current simulation time
    SimulationTime *p_simulation_time = SimulationTime::Instance();
    double time = p_simulation_time->GetDimensionalisedTime();
    
    // Set up cell type counter
    unsigned cell_type_counter[mCellTypeCount.size()];
    for (unsigned i=0; i<NUM_CELL_TYPES; i++)
    {
        cell_type_counter[i] = 0;
    }
    
    // Set up cell mutation state counter
    unsigned cell_mutation_state_counter[mCellMutationStateCount.size()];
    for (unsigned i=0; i<NUM_CELL_MUTATION_STATES; i++)
    {
        cell_mutation_state_counter[i] = 0;
    }
    
    // Set up cell cycle phase counter
    unsigned cell_cycle_phase_counter[5];
    for (unsigned i=0; i<5; i++)
    {
        cell_cycle_phase_counter[i] = 0;
    }
    
    *mpNodeFile <<  time << "\t";
    
    if (outputCellMutationStates)
    {
        *mpCellMutationStatesFile <<  time << "\t";
    }
    
    if (outputCellTypes)
    {
        *mpCellTypesFile <<  time << "\t";
    }
    
    if (outputCellVariables)
    {
        *mpCellVariablesFile <<  time << "\t";
    }
    
    if (outputCellCyclePhases)
    {
        *mpCellCyclePhasesFile <<  time << "\t";
    }

    // Write node data to file
    for (unsigned index=0; index<GetNumNodes(); index++)
    {
        unsigned colour = STEM_COLOUR; // all green if no cells have been passed in
        
        std::vector<double> proteins; // only used if outputCellVariables = true
        
        if (IsGhostNode(index) == true)
        {
            colour = INVISIBLE_COLOUR;
        }
        else if (GetNode(index)->IsDeleted())
        {
            // Do nothing
        }
        else 
        {
            TissueCell* p_cell = mNodeCellMap[index];
            
            if (outputCellCyclePhases)
            {
                switch (p_cell->GetCellCycleModel()->GetCurrentCellCyclePhase())
                {
                    case G_ZERO_PHASE:
                        cell_cycle_phase_counter[0]++;                        
                        break;                
                    case G_ONE_PHASE:
                        cell_cycle_phase_counter[1]++;                        
                        break;
                    case S_PHASE:
                        cell_cycle_phase_counter[2]++;                        
                        break;
                    case G_TWO_PHASE:
                        cell_cycle_phase_counter[3]++;                        
                        break;
                     case M_PHASE:
                        cell_cycle_phase_counter[4]++;                        
                        break;           
                    default:
                        NEVER_REACHED;
                }
            }
            
            if (mNodeCellMap[index]->GetAncestor()!=UNSIGNED_UNSET)
            {
                colour = SPECIAL_LABEL_START + p_cell->GetAncestor();
            }
            else if (mCells.size() > 0)
            {
                CellMutationState mutation = p_cell->GetMutationState();
                
                // Set colours dependent on cell type
                switch (p_cell->GetCellType())
                {
                    case STEM:
                        colour = STEM_COLOUR;
                        if (outputCellTypes)
                        {
                            cell_type_counter[0]++;
                        }
                        break;
                    case TRANSIT:
                        colour = TRANSIT_COLOUR;
                        if (outputCellTypes)
                        {
                            cell_type_counter[1]++;
                        }
                        break;
                    case DIFFERENTIATED:
                        colour = DIFFERENTIATED_COLOUR;
                        if (outputCellTypes)
                        {
                            cell_type_counter[2]++;
                        }
                        break;
                    case NECROTIC:
                        colour = APOPTOSIS_COLOUR; // paint necrotic and apoptotic cells the same colour
                        if (outputCellTypes)
                        {
                            cell_type_counter[3]++;
                        }
                        break;    
                    default:
                        NEVER_REACHED;
                }
            
                // Override colours for mutant or labelled cells
                switch (mutation)
                {
                    case HEALTHY:
                        if (outputCellMutationStates)
                        {
                            cell_mutation_state_counter[0]++;
                        }
                        break;                
                    case APC_ONE_HIT:
                        colour = EARLY_CANCER_COLOUR;
                        if (outputCellMutationStates)
                        {
                            cell_mutation_state_counter[2]++;
                        }
                        break;
                    case APC_TWO_HIT:
                        colour = LATE_CANCER_COLOUR;
                        if (outputCellMutationStates)
                        {
                            cell_mutation_state_counter[3]++;
                        }
                        break;
                    case BETA_CATENIN_ONE_HIT:
                        colour = LATE_CANCER_COLOUR;
                        if (outputCellMutationStates)
                        {
                            cell_mutation_state_counter[4]++;
                        }
                        break;
                    case LABELLED:
                        colour = LABELLED_COLOUR;
                        if (outputCellMutationStates)
                        {
                            cell_mutation_state_counter[1]++;
                        }
                        break;
                    default:
                        NEVER_REACHED;
                }
                
                if (p_cell->HasApoptosisBegun())
                {   
                    // For any type of cell set the colour to this if it is undergoing apoptosis.
                    colour = APOPTOSIS_COLOUR;   
                }
                
                if (outputCellVariables)
                {
                    proteins = p_cell->GetCellCycleModel()->GetProteinConcentrations();
                }
            }
        }
        
        if ( !(GetNode(index)->IsDeleted()) )
        {
            const c_vector<double,DIM>& position = GetNode(index)->rGetLocation();
            
            for (unsigned i=0; i<DIM; i++)
            {
                *mpNodeFile << position[i] << " ";
            }
            *mpNodeFile << colour << " ";
            
            // Write cell variable data to file if required
            if (outputCellVariables)
            {
                // Loop over cell positions
                for (unsigned i=0; i<DIM; i++)
                {
                    *mpCellVariablesFile << position[i] << " ";
                }
                // Loop over cell variables
                for (unsigned i=0; i<proteins.size(); i++)
                {
                    *mpCellVariablesFile << proteins[i] << " " ;
                }
            }
        }
    }
    
    *mpNodeFile << "\n";
   
    // Write cell mutation state data to file if required
    if (outputCellMutationStates)
    {
        for (unsigned i=0; i<NUM_CELL_MUTATION_STATES; i++)
        {
            mCellMutationStateCount[i] = cell_mutation_state_counter[i];
            *mpCellMutationStatesFile <<  cell_mutation_state_counter[i] << "\t";
        }
        *mpCellMutationStatesFile <<  "\n";
    }
    
    // Write cell type data to file if required
    if (outputCellTypes)
    {
        for (unsigned i=0; i<NUM_CELL_TYPES; i++)
        {
            mCellTypeCount[i] = cell_type_counter[i];
            *mpCellTypesFile <<  cell_type_counter[i] << "\t";
        }
        *mpCellTypesFile <<  "\n";
    }
    
    if (outputCellVariables)
    {
        *mpCellVariablesFile <<  "\n";
    }
    
    // Write cell cycle phase data to file if required
    if (outputCellCyclePhases)
    {
        for (unsigned i=0; i<5; i++)
        {
            mCellCyclePhaseCount[i] = cell_cycle_phase_counter[i];
            *mpCellCyclePhasesFile <<  cell_cycle_phase_counter[i] << "\t";            
        }
        *mpCellCyclePhasesFile <<  "\n";
    }
}

#endif //ABSTRACTTISSUE_CPP
