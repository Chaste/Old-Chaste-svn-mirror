#include "StochasticCellCycleModel.hpp"
#include "CancerParameters.hpp"


AbstractCellCycleModel *StochasticCellCycleModel::CreateCellCycleModel()
{
    return new StochasticCellCycleModel(mpGen);
}

bool StochasticCellCycleModel::ReadyToDivide(double timeSinceBirth)
{
    bool ready;
    
    
    CancerParameters *p_params = CancerParameters::Instance();
    
    switch (mCellType)
    {
        case STEM:
            ready = timeSinceBirth >= p_params->GetStemCellCycleTime();
            break;
        case TRANSIT:
            ready = timeSinceBirth >= mpGen->NormalRandomDeviate(p_params->GetTransitCellCycleTime(), 1.0);
            break;
        default:
            ready = false;
            break;
    }
    
    return ready;
}
