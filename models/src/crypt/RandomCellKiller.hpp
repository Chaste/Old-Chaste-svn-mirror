#ifndef RANDOMCELLKILLER_HPP_
#define RANDOMCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "RandomNumberGenerator.hpp"



template <unsigned SPACE_DIM>
class RandomCellKiller : public AbstractCellKiller<SPACE_DIM>
{

public:
    RandomCellKiller(Crypt<SPACE_DIM>* pCrypt)
        : AbstractCellKiller<SPACE_DIM>(pCrypt)
    {
    }
    void TestAndLabelSingleCellForApoptosis(MeinekeCryptCell& cell)
    {
        if (!cell.HasApoptosisBegun() &&
            RandomNumberGenerator::Instance()->ranf() > 0.95)
        {
            cell.StartApoptosis();
        }        
    }


    void TestAndLabelCellsForApoptosis()
    {
        for (typename Crypt<SPACE_DIM>::Iterator cell_iter = this->mpCrypt->Begin();
             cell_iter != this->mpCrypt->End();
             ++cell_iter)
        {
            TestAndLabelSingleCellForApoptosis(*cell_iter);
        }        
    }
    
    
    void RemoveDeadCells()
    {
        this->mpCrypt->RemoveDeadCells();
           
    }
};




#endif /*RANDOMCELLKILLER_HPP_*/
