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
             : mCells(rCells.begin(), rCells.end())
{
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
bool AbstractTissue<DIM>::GetIsGhostNode(unsigned index)
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
    return !(mrTissue.GetIsGhostNode(mNodeIndex) || GetNode()->IsDeleted() || (*this)->IsDead());
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



#endif //ABSTRACTTISSUE_CPP
