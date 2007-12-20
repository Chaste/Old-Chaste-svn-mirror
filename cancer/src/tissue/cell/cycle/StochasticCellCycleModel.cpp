#include "StochasticCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "Exception.hpp"
#include <cassert>
#include <iostream>


AbstractCellCycleModel *StochasticCellCycleModel::CreateDaughterCellCycleModel()
{
    return new StochasticCellCycleModel(mG1Duration, mGeneration);  // use a private constructor that doesn't reset mG1Duration.
}


void StochasticCellCycleModel::SetG1Duration()
{
    assert(mpCell!=NULL);
     
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance(); 
    
    switch (mpCell->GetCellType())
    {
        case STEM:
            mG1Duration = 1 + 4*p_gen->ranf(); // U[1,5] according to Meineke
            break;
        case TRANSIT:
            mG1Duration = 1 + 2*p_gen->ranf(); // U[1,3] according to Meineke
            break;
        case DIFFERENTIATED:
            mG1Duration = DBL_MAX;
            break;
        default:
            NEVER_REACHED;
    }
}
