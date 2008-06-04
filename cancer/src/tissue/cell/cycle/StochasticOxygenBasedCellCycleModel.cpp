/*

Copyright (C) University of Oxford, 2008

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
#include "StochasticOxygenBasedCellCycleModel.hpp"


void StochasticOxygenBasedCellCycleModel::SetG2Duration()
{
    CancerParameters* p_params = CancerParameters::Instance();
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    
    double mean = p_params->GetG2Duration();
    double standard_deviation = 1.0;
    
    mG2Duration = p_gen->NormalRandomDeviate(mean, standard_deviation);
    
    // Check that the normal random deviate has not returned a small or negative G2 duration
    if (mG2Duration < p_params->GetMinimumGapDuration())
    {
        #define COVERAGE_IGNORE
        mG2Duration = p_params->GetMinimumGapDuration();
        #undef COVERAGE_IGNORE
    }
}


void StochasticOxygenBasedCellCycleModel::InitialiseDaughterCell()
{
    AbstractSimpleCellCycleModel::InitialiseDaughterCell();
    SetG2Duration();
}


void StochasticOxygenBasedCellCycleModel::Initialise()
{
    AbstractSimpleCellCycleModel::Initialise();
    SetG2Duration();
}


void StochasticOxygenBasedCellCycleModel::ResetForDivision()
{
    AbstractSimpleCellCycleModel::ResetForDivision();
    SetG2Duration();
}


double StochasticOxygenBasedCellCycleModel::GetG2Duration()
{
    return mG2Duration;
}


StochasticOxygenBasedCellCycleModel::StochasticOxygenBasedCellCycleModel() :
    mTimeSpentInG1Phase(0.0),
    mCurrentHypoxicDuration(0.0)
{
    mCurrentHypoxiaOnsetTime = SimulationTime::Instance()->GetDimensionalisedTime();
}


double StochasticOxygenBasedCellCycleModel::GetCurrentHypoxicDuration()
{
    return mCurrentHypoxicDuration;
}
    
    
double StochasticOxygenBasedCellCycleModel::GetCurrentHypoxiaOnsetTime()
{
    return mCurrentHypoxiaOnsetTime;
}


void StochasticOxygenBasedCellCycleModel::UpdateCellCyclePhase()
{
    // mG1Duration is set when the cell cycle model is given a cell
    
    if (this->mpCell->GetCellType()!=NECROTIC)
    {
        UpdateHypoxicDuration();
        
        // Get cell's oxygen concentration
        double oxygen_concentration = CellwiseData<2>::Instance()->GetValue(mpCell,0);

        AbstractSimpleCellCycleModel::UpdateCellCyclePhase();

        if (mCurrentCellCyclePhase == G_ONE_PHASE)
        {            
            // Update G1 duration based on oxygen concentration
            double dt = SimulationTime::Instance()->GetTimeStep();
            double quiescent_concentration = CancerParameters::Instance()->GetHepaOneCellQuiescentConcentration();
           
            if (oxygen_concentration < quiescent_concentration)
            {
                mG1Duration += (1 - std::max(oxygen_concentration,0.0)/quiescent_concentration)*dt;
                mTimeSpentInG1Phase += dt;
            }
	    }
    }
}


AbstractCellCycleModel* StochasticOxygenBasedCellCycleModel::CreateDaughterCellCycleModel()
{
    return new StochasticOxygenBasedCellCycleModel(mG1Duration, mGeneration, mCurrentHypoxicDuration, mCurrentHypoxiaOnsetTime, mG2Duration);
}


void StochasticOxygenBasedCellCycleModel::UpdateHypoxicDuration()
{
    assert(this->mpCell->GetCellType()!=NECROTIC);
    assert(!this->mpCell->HasApoptosisBegun());
    
    double oxygen_concentration = CellwiseData<2>::Instance()->GetValue(this->mpCell);
    double hypoxic_concentration = CancerParameters::Instance()->GetHepaOneCellHypoxicConcentration();
    
    if ( oxygen_concentration < hypoxic_concentration)
    {
        // Update the duration of the current period of hypoxia
        mCurrentHypoxicDuration = (SimulationTime::Instance()->GetDimensionalisedTime() - mCurrentHypoxiaOnsetTime);      
                
        // Include a little bit of stochasticity here
        double prob_of_death = 0.9 - 0.5*(oxygen_concentration/hypoxic_concentration); 
        if (mCurrentHypoxicDuration > CancerParameters::Instance()->GetCriticalHypoxicDuration() && RandomNumberGenerator::Instance()->ranf() < prob_of_death)
        {                     
            this->mpCell->SetCellType(NECROTIC);
        }
    }
    else
    {
        // Reset the cell's hypoxic duration and update the time at which the onset of hypoxia occurs
        mCurrentHypoxicDuration = 0.0;
        mCurrentHypoxiaOnsetTime = SimulationTime::Instance()->GetDimensionalisedTime();
    }       
} 
