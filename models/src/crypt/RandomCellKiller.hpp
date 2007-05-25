#ifndef RANDOMCELLKILLER_HPP_
#define RANDOMCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "RandomNumberGenerator.hpp"



template <unsigned SPACE_DIM>
class RandomCellKiller : public AbstractCellKiller<SPACE_DIM>
{
private:
    double mProbability;
public:
    RandomCellKiller(Crypt<SPACE_DIM>* pCrypt, double probability)
        : AbstractCellKiller<SPACE_DIM>(pCrypt),
          mProbability(probability)
    {
        if((mProbability<0) || (mProbability>1))
        {
            EXCEPTION("Probabilitu of death must be between zero and one");
        }
    }
    
    void TestAndLabelSingleCellForApoptosis(MeinekeCryptCell& cell)
    {
        if (!cell.HasApoptosisBegun() &&
            RandomNumberGenerator::Instance()->ranf() < mProbability)
        {
            cell.StartApoptosis();
        }        
    }

    virtual void TestAndLabelCellsForApoptosis()
    {
        for (typename Crypt<SPACE_DIM>::Iterator cell_iter = this->mpCrypt->Begin();
             cell_iter != this->mpCrypt->End();
             ++cell_iter)
        {
            TestAndLabelSingleCellForApoptosis(*cell_iter);
        }        
    }
};




#endif /*RANDOMCELLKILLER_HPP_*/
