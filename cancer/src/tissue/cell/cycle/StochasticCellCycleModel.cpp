#include "StochasticCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "Exception.hpp"
#include <cassert>


AbstractCellCycleModel *StochasticCellCycleModel::CreateCellCycleModel()
{
    return new StochasticCellCycleModel(mG1Duration);  // use a private constructor that doesn't reset mG1Duration.
}

void StochasticCellCycleModel::ResetModel()
{
    mBirthTime = SimulationTime::Instance()->GetDimensionalisedTime();
    SetG1Duration();
}

void StochasticCellCycleModel::SetCell(TissueCell* pCell)
{
    AbstractCellCycleModel::SetCell(pCell);
    SetG1Duration();
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
        case HEPA_ONE:
            mG1Duration = 1 + 4*p_gen->ranf(); // U[1,5] according to Meineke
            break;
        case DIFFERENTIATED:
            mG1Duration = DBL_MAX;
            break;
        default:
            NEVER_REACHED;
    }
}

bool StochasticCellCycleModel::ReadyToDivide()
{    
    bool ready = false;
    
    CancerParameters *p_params = CancerParameters::Instance();
    double time_since_birth = GetAge();
    assert(time_since_birth>=0);
    if (mpCell->GetCellType()==DIFFERENTIATED)
    {
        mCurrentCellCyclePhase = G_ZERO;   
    }
    else if ( time_since_birth < p_params->GetMDuration() )
    {
        mCurrentCellCyclePhase = M;   
    }
    else if ( time_since_birth < p_params->GetMDuration() + mG1Duration)
    {
        mCurrentCellCyclePhase = G_ONE;   
    }
    else if ( time_since_birth < p_params->GetMDuration() + mG1Duration + p_params->GetSDuration())
    {
        mCurrentCellCyclePhase = S;   
    }
    else if ( time_since_birth < p_params->GetMDuration() + mG1Duration + p_params->GetSDuration()  + p_params->GetG2Duration())
    {
        mCurrentCellCyclePhase = G_TWO;   
    }
    else
    {
        ready = true;
        mCurrentCellCyclePhase = M;
    }
    
    return ready;
}
