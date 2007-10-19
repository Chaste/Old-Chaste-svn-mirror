#include "StochasticCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "Exception.hpp"
#include <cassert>


AbstractCellCycleModel *StochasticCellCycleModel::CreateCellCycleModel()
{
    return new StochasticCellCycleModel(mDivisionAge);  // use a private constructor that doesn't reset mDivisionAge.
}

void StochasticCellCycleModel::ResetModel()
{
    mBirthTime = SimulationTime::Instance()->GetDimensionalisedTime();
    SetDivisionAge();
}

void StochasticCellCycleModel::SetCell(TissueCell* pCell)
{
    AbstractCellCycleModel::SetCell(pCell);
    SetDivisionAge();
}

void StochasticCellCycleModel::SetDivisionAge()
{
    assert(mpCell!=NULL);
    
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    CancerParameters* p_params = CancerParameters::Instance(); 
    
    switch (mpCell->GetCellType())
    {
        case STEM:
            mDivisionAge = p_params->GetStemCellG1Duration() + p_params->GetSG2MDuration();
            break;
        case TRANSIT:
            mDivisionAge = p_gen->NormalRandomDeviate(p_params->GetTransitCellG1Duration() + p_params->GetSG2MDuration(), 1.0);
            break;
        case DIFFERENTIATED:
            mDivisionAge = DBL_MAX;
            break;
        default:
            NEVER_REACHED;
    }
}


bool StochasticCellCycleModel::ReadyToDivide()
{
    bool ready = false;
    
    if (GetAge()>=mDivisionAge)
    {
        ready = true;
    }
    
    return ready;
}
