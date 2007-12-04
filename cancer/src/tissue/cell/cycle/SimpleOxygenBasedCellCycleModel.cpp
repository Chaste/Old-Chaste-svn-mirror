#include "SimpleOxygenBasedCellCycleModel.hpp"

SimpleOxygenBasedCellCycleModel::SimpleOxygenBasedCellCycleModel() :
      mTimeSpentInG1Phase(0.0),
      mHypoxicDuration(0.0)
{
    mHypoxicDurationUpdateTime = SimulationTime::Instance()->GetDimensionalisedTime();
}


double SimpleOxygenBasedCellCycleModel::GetHypoxicDuration()
{
    return mHypoxicDuration;
}
    

bool SimpleOxygenBasedCellCycleModel::ReadyToDivide()
{
    CancerParameters *p_params = CancerParameters::Instance();
    
    // mG1Duration is set when the cell cycle model is given a cell
    
	bool ready = false;
    
    if (this->mpCell->GetCellType()!=NECROTIC)
    {
        UpdateHypoxicDuration();
        
        // get cell's oxygen concentration
        double oxygen_concentration = CellwiseData<2>::Instance()->GetValue(mpCell,0);
                                
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
            }
        }
    }    
    return ready;
}


AbstractCellCycleModel* SimpleOxygenBasedCellCycleModel::CreateCellCycleModel()
{
    return new SimpleOxygenBasedCellCycleModel(mG1Duration);
}


void SimpleOxygenBasedCellCycleModel::UpdateHypoxicDuration()
{
    assert(this->mpCell->GetCellType()!=NECROTIC);
    assert(!this->mpCell->HasApoptosisBegun());
    
    double oxygen_concentration = CellwiseData<2>::Instance()->GetValue(this->mpCell);
            
    // we want the oxygen concentration to be positive, to within numerical tolerances (hence the -1e-8)
    // assert(oxygen_concentration >= -1e-8);
            
    if ( oxygen_concentration < CancerParameters::Instance()->GetHepaOneCellHypoxicConcentration() )
    {
        // add necessary time interval to the hypoxic duration
        mHypoxicDuration += SimulationTime::Instance()->GetDimensionalisedTime() - mHypoxicDurationUpdateTime;
                
        // update mHypoxicDurationUpdateTime
        mHypoxicDurationUpdateTime = SimulationTime::Instance()->GetDimensionalisedTime();
                
        // a little bit of stochasticity here
        double prob_of_death = 0.9 - 0.5*(oxygen_concentration/CancerParameters::Instance()->GetHepaOneCellHypoxicConcentration()); 
          
        if (mHypoxicDuration > CancerParameters::Instance()->GetCriticalHypoxicDuration() && RandomNumberGenerator::Instance()->ranf() < prob_of_death)
        {                     
            this->mpCell->SetCellType(NECROTIC);
        }
    }
    else
    {
        // reset the cell's hypoxic duration
        mHypoxicDuration = 0.0;
    }       
} 
