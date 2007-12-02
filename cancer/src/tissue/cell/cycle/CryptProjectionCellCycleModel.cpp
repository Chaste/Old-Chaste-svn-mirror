#include "CryptProjectionCellCycleModel.hpp"
#include "WntGradient.hpp"
#include "Exception.hpp"
#include <iostream>
#include <cassert>

AbstractCellCycleModel *CryptProjectionCellCycleModel::CreateCellCycleModel()
{
    // use a private constructor that doesn't reset mG1Duration.
    return new CryptProjectionCellCycleModel(mG1Duration);  
}

void CryptProjectionCellCycleModel::SetG1Duration()
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
        default:
            NEVER_REACHED;
    }
}

bool CryptProjectionCellCycleModel::ReadyToDivide()
{
    assert(mpCell!=NULL);
    bool ready = false;
    
    CancerParameters *p_params = CancerParameters::Instance();
    
    double wnt_division_threshold = DBL_MAX;
    
    // set up under what level of Wnt stimulus a cell will divide
    switch (mpCell->GetMutationState())
    {
        case HEALTHY:
            wnt_division_threshold = CancerParameters::Instance()->GetRadialWntThreshold();
            break;
        case LABELLED:
            wnt_division_threshold = CancerParameters::Instance()->GetRadialWntThreshold();
            break;
        default:
            NEVER_REACHED;
    }
    
    // get the age of the cell
    double time_since_birth = GetAge();
    assert(time_since_birth>=0);
        
    // update the current cell cycle phase    
    if ( time_since_birth < p_params->GetMDuration() )
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
    	if (WntGradient::Instance()->GetWntLevel(mpCell) > wnt_division_threshold)
        {
            mpCell->SetCellType(STEM);
        }
        else
        {
            mpCell->SetCellType(TRANSIT);
        }
	    ready = true;
	    mCurrentCellCyclePhase = M;
	}    
    
    return ready;
}
