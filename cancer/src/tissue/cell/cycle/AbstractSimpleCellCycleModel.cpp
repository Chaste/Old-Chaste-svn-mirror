#include "AbstractSimpleCellCycleModel.hpp"

double AbstractSimpleCellCycleModel::GetG1Duration()
{
	return mG1Duration;	
}

void AbstractSimpleCellCycleModel::SetG1Duration()
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


void AbstractSimpleCellCycleModel::SetCell(TissueCell* pCell)
{
    AbstractCellCycleModel::SetCell(pCell);
    // This method should only be called once per cell cycle model - when it is created so G1Duration can be set here.
    SetG1Duration();	
}

void AbstractSimpleCellCycleModel::ResetModel()
{
    mBirthTime = SimulationTime::Instance()->GetDimensionalisedTime();
    SetG1Duration();
}


bool AbstractSimpleCellCycleModel::ReadyToDivide()
{
    assert(mpCell != NULL);
    bool ready = false;
    
    double time_since_birth = GetAge();
    assert(time_since_birth>=0);
    
    if (mpCell->GetCellType()==DIFFERENTIATED)
    {
        mCurrentCellCyclePhase = G_ZERO_PHASE;   
    }
    else if ( time_since_birth < GetMDuration() )
    {
        mCurrentCellCyclePhase = M_PHASE;   
    }
    else if ( time_since_birth < GetMDuration() + mG1Duration)
    {
        mCurrentCellCyclePhase = G_ONE_PHASE;   
    }
    else if ( time_since_birth < GetMDuration() + mG1Duration + GetSDuration())
    {
        mCurrentCellCyclePhase = S_PHASE;   
    }
    else if ( time_since_birth < GetMDuration() + mG1Duration + GetSDuration()  + GetG2Duration())
    {
        mCurrentCellCyclePhase = G_TWO_PHASE;   
    }
    else
    {
        ready = true;
        mCurrentCellCyclePhase = M_PHASE;
    }
    
    return ready;
}


