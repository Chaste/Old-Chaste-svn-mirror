#ifndef RANDOMCELLKILLER_HPP_
#define RANDOMCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "RandomNumberGenerator.hpp"


/**
 *  Randomly kills cells based on the user set probability
 *  The probability passed into the constructor will be the probability
 *  of any cell dying whenever this TestAndLabelCellsForApoptosis is called.
 *  Note this does take into account current times or timesteps, so if
 *  more timesteps are used, and TestAndLabelCellsForApoptosis() is callled 
 *  at each timestep, more cells will die.  
 */
template <unsigned SPACE_DIM>
class RandomCellKiller : public AbstractCellKiller<SPACE_DIM>
{
private:
    double mProbabilityOfDeath;
public:
    RandomCellKiller(Crypt<SPACE_DIM>* pCrypt, double probabilityOfDeath)
        : AbstractCellKiller<SPACE_DIM>(pCrypt),
          mProbabilityOfDeath(probabilityOfDeath)
    {
        if((mProbabilityOfDeath<0) || (mProbabilityOfDeath>1))
        {
            EXCEPTION("Probabilitu of death must be between zero and one");
        }
    }
    
    void TestAndLabelSingleCellForApoptosis(MeinekeCryptCell& cell)
    {
        if (!cell.HasApoptosisBegun() &&
            RandomNumberGenerator::Instance()->ranf() < mProbabilityOfDeath)
        {
            cell.StartApoptosis();
        }        
    }

    /**
     *  Loops over cells and starts apoptosis randomly, based on the user-set 
     *  probability
     */
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
