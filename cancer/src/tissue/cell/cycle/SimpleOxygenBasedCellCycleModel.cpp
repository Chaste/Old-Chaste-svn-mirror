#include "SimpleOxygenBasedCellCycleModel.hpp"

SimpleOxygenBasedCellCycleModel::SimpleOxygenBasedCellCycleModel() :
      mTimeSpentInG1Phase(0.0),
      mCurrentHypoxicDuration(0.0)
{
    mCurrentHypoxiaOnsetTime = SimulationTime::Instance()->GetDimensionalisedTime();
}


double SimpleOxygenBasedCellCycleModel::GetCurrentHypoxicDuration()
{
    return mCurrentHypoxicDuration;
}
    
    
double SimpleOxygenBasedCellCycleModel::GetCurrentHypoxiaOnsetTime()
{
    return mCurrentHypoxiaOnsetTime;
}


bool SimpleOxygenBasedCellCycleModel::ReadyToDivide()
{
    // mG1Duration is set when the cell cycle model is given a cell
    
    bool ready = false;
	  
    if (this->mpCell->GetCellType()!=NECROTIC)
    {
        UpdateHypoxicDuration();
        
        // get cell's oxygen concentration
        double oxygen_concentration = CellwiseData<2>::Instance()->GetValue(mpCell,0);

	    ready = AbstractSimpleCellCycleModel::ReadyToDivide();

        if (mCurrentCellCyclePhase == G_ONE_PHASE)
        {
	       // Update G1 duration based on oxygen concentration
	       double dt = SimulationTime::Instance()->GetTimeStep();
	       mG1Duration += (1-std::max(oxygen_concentration,0.0))*dt;
	       mTimeSpentInG1Phase += dt;
	    }
    }    
    return ready;
}


AbstractCellCycleModel* SimpleOxygenBasedCellCycleModel::CreateDaughterCellCycleModel()
{
    return new SimpleOxygenBasedCellCycleModel(mG1Duration, mGeneration, mCurrentHypoxicDuration, mCurrentHypoxiaOnsetTime);
}


void SimpleOxygenBasedCellCycleModel::UpdateHypoxicDuration()
{
    assert(this->mpCell->GetCellType()!=NECROTIC);
    assert(!this->mpCell->HasApoptosisBegun());
    
    double oxygen_concentration = CellwiseData<2>::Instance()->GetValue(this->mpCell);
            
    if ( oxygen_concentration < CancerParameters::Instance()->GetHepaOneCellHypoxicConcentration() )
    {
        // update the duration of the current period of hypoxia
        mCurrentHypoxicDuration = (SimulationTime::Instance()->GetDimensionalisedTime() - mCurrentHypoxiaOnsetTime);      
                
        // a little bit of stochasticity here
        double prob_of_death = 0.9 - 0.5*(oxygen_concentration/CancerParameters::Instance()->GetHepaOneCellHypoxicConcentration()); 
        if (mCurrentHypoxicDuration > CancerParameters::Instance()->GetCriticalHypoxicDuration() && RandomNumberGenerator::Instance()->ranf() < prob_of_death)
        {                     
            this->mpCell->SetCellType(NECROTIC);
        }
    }
    else
    {
        // reset the cell's hypoxic duration and update the time at which the onset of hypoxia occurs
        mCurrentHypoxicDuration = 0.0;
        mCurrentHypoxiaOnsetTime = SimulationTime::Instance()->GetDimensionalisedTime();
    }       
} 
