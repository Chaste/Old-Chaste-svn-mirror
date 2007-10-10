#include "StochasticCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "Exception.hpp"
#include <cassert>


AbstractCellCycleModel *StochasticCellCycleModel::CreateCellCycleModel()
{
    return new StochasticCellCycleModel();
}


void StochasticCellCycleModel::ResetModel()
{
    mBirthTime = SimulationTime::Instance()->GetDimensionalisedTime();
}

bool StochasticCellCycleModel::ReadyToDivide()
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    CancerParameters* p_params = CancerParameters::Instance(); 
    
    bool ready;
    
    double timeSinceBirth = GetAge();
    
    switch (mpCell->GetCellType())
    {
        case STEM:
            ready = (timeSinceBirth >= p_params->GetStemCellCycleTime());
            break;
        case TRANSIT:
            ready = (timeSinceBirth >= p_gen->NormalRandomDeviate(p_params->GetTransitCellCycleTime(), 1.0));
            break;
        case HEPA_ONE:
            ready = (timeSinceBirth >= p_gen->NormalRandomDeviate(p_params->GetHepaOneCellCycleTime(), 1.0));
            break;
        default:
            ready = false;
            break;
    }
    
    return ready;
}
