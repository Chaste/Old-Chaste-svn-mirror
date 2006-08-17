#include "FixedCellCycleModel.hpp"
#include "CancerParameters.hpp"

AbstractCellCycleModel *FixedCellCycleModel::CreateCellCycleModel()
{
    return new FixedCellCycleModel();
}

bool FixedCellCycleModel::ReadyToDivide(double timeSinceBirth)
{
    bool ready;
    
    CancerParameters *p_params = CancerParameters::Instance();
    
    switch (mCellType)
    {
        case STEM:
            ready = timeSinceBirth >= p_params->GetStemCellCycleTime();
            break;
        case TRANSIT:
            ready = timeSinceBirth >= p_params->GetTransitCellCycleTime();
            break;
        default:
            ready = false;
            break;
    }
    
    return ready;
}
