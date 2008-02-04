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
TissueCell& AbstractTissue<DIM>::rGetCellAtNodeIndex(unsigned index)
{
    return *(mNodeCellMap[index]);
}


#endif //ABSTRACTTISSUE_CPP
