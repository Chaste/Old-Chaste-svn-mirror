#include "StochasticCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "Exception.hpp"
#include <cassert>


AbstractCellCycleModel *StochasticCellCycleModel::CreateCellCycleModel()
{
    return new StochasticCellCycleModel();
}


StochasticCellCycleModel::StochasticCellCycleModel()
{
    SimulationTime* p_sim_time = SimulationTime::Instance();
    if (p_sim_time->IsStartTimeSetUp()==false)
    {
        EXCEPTION("StochasticCellCycleModel is being created but SimulationTime has not been set up");
    }  
    mBirthTime = p_sim_time->GetDimensionalisedTime();
}

void StochasticCellCycleModel::SetBirthTime(double birthTime)
{
    mBirthTime = birthTime;
}

void StochasticCellCycleModel::ResetModel()
{
    mBirthTime = SimulationTime::Instance()->GetDimensionalisedTime();
}

bool StochasticCellCycleModel::ReadyToDivide(std::vector<double> cellCycleInfluences)
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    CancerParameters* p_params = CancerParameters::Instance(); 
    //assert(cellCycleInfluences.size()==0);
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
        default:
            ready = false;
            break;
    }
    
    return ready;
}
