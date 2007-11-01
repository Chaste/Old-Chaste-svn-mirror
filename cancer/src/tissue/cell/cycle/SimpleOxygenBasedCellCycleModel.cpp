#include "SimpleOxygenBasedCellCycleModel.hpp"

bool SimpleOxygenBasedCellCycleModel::ReadyToDivide()
{
    CancerParameters *p_params = CancerParameters::Instance();
    
    // mG1Duration is now set when the cell cycle model is given a cell.
    
	bool ready = false;
    
    // get cell's oxygen concentration
    double oxygen_concentration = CellwiseData<2>::Instance()->GetValue(mpCell);
    
    // we want the oxygen concentration to be positive,
    // to within numerical tolerances (hence the -1e-8)
//        assert(oxygen_concentration >= -1e-8);
                    
    double time_since_birth = GetAge();
           
    if (mpCell->GetCellType()==DIFFERENTIATED)
    {
        mCurrentCellCyclePhase = G_ZERO;   
    }
    else 
    {
        if ( GetAge() < p_params->GetMDuration() )
        {
            mCurrentCellCyclePhase = M;   
        }
        else if ( time_since_birth < p_params->GetMDuration() + mG1Duration )
        {
            mCurrentCellCyclePhase = G_ONE;
             
            mG1Duration = mG1Duration + (1-std::max(oxygen_concentration,0.0))*SimulationTime::Instance()->GetTimeStep();
            mTimeSpentInG1Phase = mTimeSpentInG1Phase + SimulationTime::Instance()->GetTimeStep();  
        }
        else if ( time_since_birth < p_params->GetMDuration() + mG1Duration + p_params->GetSDuration() )
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
            // mCurrentCellCyclePhase = M;
        }
    }
           
    return ready;
    
}

AbstractCellCycleModel* SimpleOxygenBasedCellCycleModel::CreateCellCycleModel()
{
    return new SimpleOxygenBasedCellCycleModel(mG1Duration);
}


