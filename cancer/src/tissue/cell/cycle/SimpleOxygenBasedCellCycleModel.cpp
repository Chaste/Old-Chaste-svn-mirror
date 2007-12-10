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
    // mG1Duration is set when the cell cycle model is given a cell
    
    bool ready = false;
	  
    if (this->mpCell->GetMutationState()!=NECROTIC)
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


AbstractCellCycleModel* SimpleOxygenBasedCellCycleModel::CreateCellCycleModel()
{
    return new SimpleOxygenBasedCellCycleModel(mG1Duration, mGeneration, mHypoxicDuration, mHypoxicDurationUpdateTime);
}


void SimpleOxygenBasedCellCycleModel::UpdateHypoxicDuration()
{
    assert(this->mpCell->GetMutationState()!=NECROTIC);
    assert(!this->mpCell->HasApoptosisBegun());
    
    double oxygen_concentration = CellwiseData<2>::Instance()->GetValue(this->mpCell);
            
    // we want the oxygen concentration to be positive, to within numerical tolerances (hence the -1e-8)
    // assert(oxygen_concentration >= -1e-8);
            
    if ( oxygen_concentration < CancerParameters::Instance()->GetHepaOneCellHypoxicConcentration() )
    {
        // add necessary time interval to the hypoxic duration
       
        mHypoxicDuration = (SimulationTime::Instance()->GetDimensionalisedTime() - mHypoxicDurationUpdateTime);      
        //Not this mHypoxicDurationUpdateTime is currently redundant.  This change has been made to avoid
        //floating point inaccuracies growing in the computation of mHypoxicDuration
        
        // DO NOT update mHypoxicDurationUpdateTime
        //mHypoxicDurationUpdateTime = SimulationTime::Instance()->GetDimensionalisedTime();
                
        // a little bit of stochasticity here
        double prob_of_death = 0.9 - 0.5*(oxygen_concentration/CancerParameters::Instance()->GetHepaOneCellHypoxicConcentration()); 
        if (mHypoxicDuration > CancerParameters::Instance()->GetCriticalHypoxicDuration() && RandomNumberGenerator::Instance()->ranf() < prob_of_death)
        {                     
            this->mpCell->SetMutationState(NECROTIC);
        }
    }
    else
    {
        // reset the cell's hypoxic duration
        mHypoxicDuration = 0.0;
    }       
} 
