#include "StochasticDivisionRuleCellCycleModel.hpp"

void StochasticDivisionRuleCellCycleModel::SetG1Duration()    
{
    assert(mpCell!=NULL);
    
    CancerParameters* p_params = CancerParameters::Instance(); 
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance(); 
    
    switch (mpCell->GetCellType())
    {
        case STEM:
            mG1Duration = p_gen->NormalRandomDeviate(p_params->GetStemCellG1Duration(),1.0);            
            break;
        case TRANSIT:
            mG1Duration = p_gen->NormalRandomDeviate(p_params->GetTransitCellG1Duration(),1.0);
            break;
        case DIFFERENTIATED:
            mG1Duration = DBL_MAX;
            break;
        default:
            NEVER_REACHED;
    }        
}

void StochasticDivisionRuleCellCycleModel::ResetForDivision()
{
    if (mGeneration+1u > CancerParameters::Instance()->GetMaxTransitGenerations())
    {
        mpCell->SetCellType(DIFFERENTIATED);
    }
    
    // If dealing with a stem cell, we may have symmetric division.
    // We therefore neglect the possibility of de-differentiation.
    // NB. This code must be implemented before the call to 
    // AbstractSimpleCellCycleModel::ResetForDivision(), because that 
    // method sets the G1 duration based on the cell type.
    if (mpCell->GetCellType() == STEM)
    {
        double test_number = RandomNumberGenerator::Instance()->ranf(); // U(0,1)
        double sym_div_prob = CancerParameters::Instance()->GetSymmetricDivisionProbability();
                    
        // If undergoing symmetric division...
        if (test_number < sym_div_prob)
        {
            mDividedSymmetrically = true;
            
            // Check if the daughter cells are both STEM or TRANSIT.
            // We assign an equal probability to each of these events.                
            if (test_number < 0.5*sym_div_prob)
            {                    
                mpCell->SetCellType(STEM);
            }
            else
            {                    
                mpCell->SetCellType(TRANSIT);
            }
        }
        else
        {
            mDividedSymmetrically = false;
        }
    }
            
    AbstractSimpleCellCycleModel::ResetForDivision();
    
    if (mpCell->GetCellType() == STEM)
    {
        mGeneration = 0;
    }
}

void StochasticDivisionRuleCellCycleModel::InitialiseDaughterCell()
{
    // If the cell was born out of symmetric division, 
    // then do not alter generation or cell type        
    if (mDividedSymmetrically == false)
    {
        if (mGeneration == 0)
        {
            mGeneration = 1;
        }
        // Daughter cell is always a TRANSIT or DIFFERENTIATED
        mpCell->SetCellType(TRANSIT);
        if (mGeneration > CancerParameters::Instance()->GetMaxTransitGenerations())
        {
            mpCell->SetCellType(DIFFERENTIATED);
        }
    }
    AbstractSimpleCellCycleModel::InitialiseDaughterCell();
}

AbstractCellCycleModel* StochasticDivisionRuleCellCycleModel::CreateDaughterCellCycleModel()
{
    // Use a private constructor that doesn't reset mG1Duration
    return new StochasticDivisionRuleCellCycleModel(mG1Duration, mGeneration, mDividedSymmetrically);
}

bool StochasticDivisionRuleCellCycleModel::DividedSymmetrically()
{
    return mDividedSymmetrically;
}
