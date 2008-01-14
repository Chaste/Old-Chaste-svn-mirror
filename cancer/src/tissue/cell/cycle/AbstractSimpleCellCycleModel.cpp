#include "AbstractSimpleCellCycleModel.hpp"

void AbstractSimpleCellCycleModel::InitialiseDaughterCell()
{
    AbstractCellCycleModel::InitialiseDaughterCell();
    SetG1Duration();
}


void AbstractSimpleCellCycleModel::Initialise() 
{ 
    SetG1Duration(); 
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
        case DIFFERENTIATED:
            mG1Duration = DBL_MAX;
            break;
        case NECROTIC:
            mG1Duration = DBL_MAX;
            break;    
        default:
            NEVER_REACHED;
    }
}


void AbstractSimpleCellCycleModel::ResetModel()
{
    AbstractCellCycleModel::ResetModel();
    mBirthTime = SimulationTime::Instance()->GetDimensionalisedTime();
    SetG1Duration();
}


void AbstractSimpleCellCycleModel::UpdateCellCyclePhase()
{
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
    else if ( time_since_birth < GetMDuration() + mG1Duration + GetSDuration() + GetG2Duration())
    {
        mCurrentCellCyclePhase = G_TWO_PHASE;   
    }    
}

