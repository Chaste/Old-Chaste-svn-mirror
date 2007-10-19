#include "SimpleWntCellCycleModel.hpp"
#include "WntGradient.hpp"
#include "Exception.hpp"
#include <iostream>
#include <cassert>

AbstractCellCycleModel *SimpleWntCellCycleModel::CreateCellCycleModel()
{
    return new SimpleWntCellCycleModel();
}

void SimpleWntCellCycleModel::ResetModel()
{
    mBirthTime = SimulationTime::Instance()->GetDimensionalisedTime();
    mCycleTime = RandomNumberGenerator::Instance()->
                        NormalRandomDeviate(
                            CancerParameters::Instance()->GetTransitCellG1Duration()
                                +  CancerParameters::Instance()->GetSG2MDuration() , 
                            1.0);   // (mean and s.d.)
}

bool SimpleWntCellCycleModel::ReadyToDivide()
{
    assert(mpCell!=NULL);
    double wnt_division_threshold = DBL_MAX;
    
    // Set up under what level of Wnt stimulus a cell will divide
    switch (mpCell->GetMutationState())
    {
        case HEALTHY:
            wnt_division_threshold = 0.5;
            break;
        case LABELLED:
            wnt_division_threshold = 0.5;
            break;
        case APC_ONE_HIT:
            wnt_division_threshold = 0.4;
            break;
        case BETA_CATENIN_ONE_HIT:
            wnt_division_threshold = 0.1;
            break;
        case APC_TWO_HIT:
            wnt_division_threshold = 0.0;
            break;
        default:
            NEVER_REACHED;
    }
    
    bool ready = false;
    CellType cell_type = DIFFERENTIATED;
    if (WntGradient::Instance()->GetWntLevel(mpCell) >= wnt_division_threshold)
    {
        cell_type = TRANSIT;
        if (GetAge() >= mCycleTime)
        {
            ready = true;   
        }
    }
    
    mpCell->SetCellType(cell_type); // update the cell type to reflect wnt conc.
    return ready;
}
