#include "FixedCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "Exception.hpp"
#include <iostream>
#include <cassert>

AbstractCellCycleModel *FixedCellCycleModel::CreateCellCycleModel()
{
    return new FixedCellCycleModel();
}

FixedCellCycleModel::FixedCellCycleModel()
{
    SimulationTime* p_sim_time = SimulationTime::Instance();
    if (p_sim_time->IsStartTimeSetUp()==false)
    {
        EXCEPTION("FixedCellCycleModel is being created but SimulationTime has not been set up");
    }
    mBirthTime = p_sim_time->GetDimensionalisedTime();
}

void FixedCellCycleModel::SetBirthTime(double birthTime)
{
    mBirthTime = birthTime;
}

void FixedCellCycleModel::ResetModel()
{
    mBirthTime = SimulationTime::Instance()->GetDimensionalisedTime();
}

bool FixedCellCycleModel::ReadyToDivide(std::vector<double> cellCycleInfluences)
{
    //assert(cellCycleInfluences.size()==0); NOT Needed - we just ignore them
    bool ready;
    
    CancerParameters *p_params = CancerParameters::Instance();
    double timeSinceBirth = GetAge();
    
    switch (mpCell->GetCellType())
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
