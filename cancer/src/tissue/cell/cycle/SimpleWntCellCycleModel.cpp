#include "SimpleWntCellCycleModel.hpp"
#include "WntGradient.hpp"
#include "Exception.hpp"
#include <iostream>
#include <cassert>

AbstractCellCycleModel *SimpleWntCellCycleModel::CreateCellCycleModel()
{
    // use a private constructor that doesn't reset mG1Duration.
    return new SimpleWntCellCycleModel(mG1Duration,mGeneration);  
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
        case HEPA_ONE:
            mG1Duration = p_gen->NormalRandomDeviate(p_params->GetHepaOneCellG1Duration(),1.0);
            break;
        case DIFFERENTIATED:
            mG1Duration = DBL_MAX;
            break;
        default:
            NEVER_REACHED;
    }        
}

bool SimpleWntCellCycleModel::ReadyToDivide()
{
    assert(mpCell!=NULL);
    bool ready = false;
    
    CancerParameters *p_params = CancerParameters::Instance();
    WntGradient* p_wnt_gradient = WntGradient::Instance();
    
    double wnt_stem_cell_threshold = DBL_MAX;
    double wnt_division_threshold = DBL_MAX;
    double healthy_threshold = 0.65;    // cell will divide if Wnt level >= to this value.
        
    // In the case of a RADIAL Wnt gradient, set up under what level 
    // of Wnt stimulus a cell will change type 
    if (p_wnt_gradient->GetType()==RADIAL)
    {
        wnt_stem_cell_threshold = p_params->GetRadialWntThreshold();
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
            wnt_division_threshold = 0.5;
            break;
        case BETA_CATENIN_ONE_HIT:  // less than above value
            wnt_division_threshold = 0.1;
            break;
        case APC_TWO_HIT:   // should be zero (no Wnt-dependence).
            wnt_division_threshold = 0.0;
            break;
        default:
            NEVER_REACHED;
    }
        
    double time_since_birth = GetAge();
    assert(time_since_birth >= 0);
    
     
    // If the Wnt stimulus is below the threshold, the cell is
    // of type DIFFERENTIATED and hence in G0 phase       
    CellType cell_type = DIFFERENTIATED;    
    mCurrentCellCyclePhase = G_ZERO_PHASE;
    
    
    // If the Wnt stimulus exceeds the threshold, the cell is
    // of type TRANSIT, and hence its cell cycle phase depends
    // on its age     
    if (p_wnt_gradient->GetWntLevel(mpCell) >= wnt_division_threshold)
    {       
        cell_type = TRANSIT;
        
        if (p_wnt_gradient->GetType()==RADIAL)
        {
            if (p_wnt_gradient->GetWntLevel(mpCell) > wnt_stem_cell_threshold)
            {
                cell_type = STEM;
            }
        }
        
        if ( time_since_birth < p_params->GetMDuration() )
        {
            mCurrentCellCyclePhase = M_PHASE;   
        }
        else if ( time_since_birth < p_params->GetMDuration() + mG1Duration)
        {
            mCurrentCellCyclePhase = G_ONE_PHASE;   
        }
        else if ( time_since_birth < p_params->GetMDuration() + mG1Duration + p_params->GetSDuration())
        {
            mCurrentCellCyclePhase = S_PHASE;   
        }
        else if ( time_since_birth < p_params->GetMDuration() + mG1Duration + p_params->GetSDuration()  + p_params->GetG2Duration())
        {
            mCurrentCellCyclePhase = G_TWO_PHASE;   
        }
        else
        {
            ready = true;          
            mCurrentCellCyclePhase = M_PHASE;
        }
    }
    
    // Update the cell type to reflect the Wnt concentration
    mpCell->SetCellType(cell_type);     
    
    return ready;
}

std::vector<CellType> SimpleWntCellCycleModel::GetNewCellTypes(CellType cellType)
{   
    std::vector<CellType> cell_types(2);
    
    if (WntGradient::Instance()->GetType()==RADIAL)
    {        
        cell_types[0] = cellType;
        cell_types[1] = TRANSIT;
        if (cellType == STEM)
        {
            SetGeneration(0u); 
        }
    }
    else
    {
        cell_types = AbstractCellCycleModel::GetNewCellTypes(cellType);
    }
    return cell_types;
}
