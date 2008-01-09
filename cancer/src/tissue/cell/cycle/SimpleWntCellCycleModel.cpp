#include "SimpleWntCellCycleModel.hpp"
#include "WntGradient.hpp"
#include "Exception.hpp"
#include <iostream>
#include <cassert>

AbstractCellCycleModel *SimpleWntCellCycleModel::CreateDaughterCellCycleModel()
{
    // use a private constructor that doesn't reset mG1Duration.
    return new SimpleWntCellCycleModel(mG1Duration, mGeneration);  
}

/** 
 * The G1 duration is taken from a normal distribution, whose mean is
 * the G1 duration given in CancerParameters for the cell type, and 
 * whose standard deviataion is 1. 
 */  
void SimpleWntCellCycleModel::SetG1Duration()
{
    assert(mpCell!=NULL);
    
    CancerParameters* p_params = CancerParameters::Instance(); 
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance(); 
    
    switch (mpCell->GetCellType())
    {
        case STEM:  // STEM cells should behave just like transit cells in a Wnt simulation
            mG1Duration = p_gen->NormalRandomDeviate(p_params->GetTransitCellG1Duration(),1.0);            
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

void SimpleWntCellCycleModel::UpdateCellCyclePhase()
{
    CancerParameters *p_params = CancerParameters::Instance();
    WntGradient* p_wnt_gradient = WntGradient::Instance();
    
    double wnt_stem_cell_threshold = DBL_MAX;
    double wnt_division_threshold = DBL_MAX;
    double healthy_threshold = p_params->GetWntTransitThreshold();    // cell will divide if Wnt level >= to this value.
        
    // In the case of a RADIAL Wnt gradient, set up under what level 
    // of Wnt stimulus a cell will change type 
    if (p_wnt_gradient->GetType()==RADIAL)
    {
        wnt_stem_cell_threshold = p_params->GetWntStemThreshold();
    }
    
    // Set up under what level of Wnt stimulus a cell will divide
    switch (mpCell->GetMutationState())
    {
        case HEALTHY:
            wnt_division_threshold = healthy_threshold;
            break;
        case LABELLED:
            wnt_division_threshold = healthy_threshold;
            break;
        case APC_ONE_HIT:   // should be less than healthy values
            wnt_division_threshold = 0.77*healthy_threshold;
            break;
        case BETA_CATENIN_ONE_HIT:  // less than above value
            wnt_division_threshold = 0.155*healthy_threshold;
            break;
        case APC_TWO_HIT:   // should be zero (no Wnt-dependence).
            wnt_division_threshold = 0.0;
            break;
        default:
            NEVER_REACHED;
    }
    
    // If the Wnt stimulus exceeds the threshold, the cell is
    // of type TRANSIT, and hence its cell cycle phase depends
    // on its age, just as in AbstractSimpleCellCycleModel.
    if (p_wnt_gradient->GetWntLevel(mpCell) >= wnt_division_threshold)
    {       
        CellType cell_type = TRANSIT;
        
        if (p_wnt_gradient->GetType()==RADIAL)
        {
            if (p_wnt_gradient->GetWntLevel(mpCell) > wnt_stem_cell_threshold)
            {
                cell_type = STEM;
            }
        }
    	// Update the cell type to reflect the Wnt concentration
    	mpCell->SetCellType(cell_type);
       
    	AbstractSimpleCellCycleModel::UpdateCellCyclePhase();
    }
    else
    {
    	// If the Wnt stimulus is below the threshold, the cell is
    	// of type DIFFERENTIATED and hence in G0 phase
    	mpCell->SetCellType(DIFFERENTIATED);
    	mCurrentCellCyclePhase = G_ZERO_PHASE;
    }
}

std::vector<CellType> SimpleWntCellCycleModel::GetNewCellTypes()
{   
    CellType cell_type = mpCell->GetCellType();
    std::vector<CellType> new_cell_types(2);
    
    if (WntGradient::Instance()->GetType()==RADIAL)
    {        
        new_cell_types[0] = cell_type;
        new_cell_types[1] = TRANSIT;
        if (cell_type == STEM)
        {
            SetGeneration(0u); 
        }
    }
    else
    {
        new_cell_types = AbstractCellCycleModel::GetNewCellTypes();
    }
    return new_cell_types;
}
