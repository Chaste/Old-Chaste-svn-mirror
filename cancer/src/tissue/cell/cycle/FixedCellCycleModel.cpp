#include "FixedCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "Exception.hpp"
#include <iostream>
#include <cassert>

AbstractCellCycleModel *FixedCellCycleModel::CreateCellCycleModel()
{
    return new FixedCellCycleModel();
}

void FixedCellCycleModel::ResetModel()
{
    mBirthTime = SimulationTime::Instance()->GetDimensionalisedTime();
}

bool FixedCellCycleModel::ReadyToDivide()
{
    bool ready = false;
    
    CancerParameters *p_params = CancerParameters::Instance();
    double timeSinceBirth = GetAge();
    
    switch (mpCell->GetCellType())
    {
        case STEM:
            ready = timeSinceBirth >= p_params->GetStemCellG1Duration() + p_params->GetSG2MDuration();
            break;
        case TRANSIT:
            ready = timeSinceBirth >= p_params->GetTransitCellG1Duration() + p_params->GetSG2MDuration();
            break;
        case HEPA_ONE:
            ready = timeSinceBirth >= p_params->GetHepaOneCellG1Duration() + p_params->GetSG2MDuration();
            break;
        default:
            ready = false;
            break;
    }
    
    return ready;
}
