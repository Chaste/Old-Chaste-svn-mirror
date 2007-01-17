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
	mpSimulationTime = SimulationTime::Instance();
	if(mpSimulationTime->IsSimulationTimeSetUp()==false)
	{
		EXCEPTION("FixedCellCycleModel is being created but SimulationTime has not been set up");
	}
	mBirthTime = mpSimulationTime->GetDimensionalisedTime();
}

void FixedCellCycleModel::ResetModel()
{
	mpSimulationTime = SimulationTime::Instance();
	mBirthTime = mpSimulationTime->GetDimensionalisedTime();	
}

bool FixedCellCycleModel::ReadyToDivide(std::vector<double> cellCycleInfluences)
{	
	mpSimulationTime = SimulationTime::Instance();
	assert(cellCycleInfluences.size()==0);
    bool ready;
    
    CancerParameters *p_params = CancerParameters::Instance();
    double timeSinceBirth = GetAge();
    
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
