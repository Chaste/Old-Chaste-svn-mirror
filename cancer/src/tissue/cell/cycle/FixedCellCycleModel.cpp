#include "FixedCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "Exception.hpp"
#include <iostream>
#include <cassert>

AbstractCellCycleModel *FixedCellCycleModel::CreateCellCycleModel()
{
    return new FixedCellCycleModel(mG1Duration);
}

void FixedCellCycleModel::ResetModel()
{
    mBirthTime = SimulationTime::Instance()->GetDimensionalisedTime();
    SetG1Duration();
}

void FixedCellCycleModel::SetCell(TissueCell* pCell)
{
    AbstractCellCycleModel::SetCell(pCell);
    SetG1Duration();
}

void FixedCellCycleModel::SetG1Duration()
{
    assert(mpCell!=NULL);
    
    CancerParameters* p_params = CancerParameters::Instance(); 
    
    switch (mpCell->GetCellType())
    {
        case STEM:
            mG1Duration = p_params->GetStemCellG1Duration();
            break;
        case TRANSIT:
            mG1Duration = p_params->GetTransitCellG1Duration();
            break;
        case HEPA_ONE:
            mG1Duration = p_params->GetHepaOneCellG1Duration();
            break;
        case DIFFERENTIATED:
            mG1Duration = DBL_MAX;
            break;
        default:
            NEVER_REACHED;
    }
}

bool FixedCellCycleModel::ReadyToDivide()
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
