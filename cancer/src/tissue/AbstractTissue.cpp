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
               mTissueContainsMesh(false)
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
    for(unsigned i=0; i<mCellTypeCount.size(); i++)
    {
        mCellTypeCount[i] = 0;
    }
}

template<unsigned DIM>
void AbstractTissue<DIM>::InitialiseCells()
{
    for(std::list<TissueCell>::iterator iter = mCells.begin();
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
unsigned AbstractTissue<DIM>::GetNumRealCells()
{
    unsigned counter = 0;
    for(typename AbstractTissue<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        counter++;
    }
    return counter;
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
c_vector<unsigned,5> AbstractTissue<DIM>::GetCellTypeCount()
{
    return mCellTypeCount;
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
    assert(mrTissue.GetGhostNodesSize() == mrTissue.GetNumNodes() );
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
void AbstractTissue<DIM>::CreateOutputFiles(const std::string &rDirectory, bool rCleanOutputDirectory, bool outputCellTypes)
{
    OutputFileHandler output_file_handler(rDirectory, rCleanOutputDirectory);
    mpNodeFile = output_file_handler.OpenOutputFile("results.viznodes");
    mpCellTypesFile = output_file_handler.OpenOutputFile("celltypes.dat");
    mpCellVariablesFile = output_file_handler.OpenOutputFile("cellvariables.dat");
    
    if (outputCellTypes)
    {
        *mpCellTypesFile <<   "Time\t Healthy\t Labelled\t APC_1\t APC_2\t BETA_CAT \n";
    }
}

template<unsigned DIM>
void AbstractTissue<DIM>::CloseOutputFiles()
{
    mpNodeFile->close();
    mpCellTypesFile->close();
    mpCellVariablesFile->close();
}

#endif //ABSTRACTTISSUE_CPP
