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
	mpSimulationTime = SimulationTime::Instance();
	if(mpSimulationTime->IsStartTimeSetUp()==false)
	{
		EXCEPTION("StochasticCellCycleModel is being created but SimulationTime has not been set up");
	}
	mpCancerParams = CancerParameters::Instance();
	mBirthTime = mpSimulationTime->GetDimensionalisedTime();

}

void StochasticCellCycleModel::SetBirthTime(double birthTime)
{
	mBirthTime = birthTime;	
}

void StochasticCellCycleModel::ResetModel()
{
	mpSimulationTime = SimulationTime::Instance();
	mBirthTime = mpSimulationTime->GetDimensionalisedTime();	
}

bool StochasticCellCycleModel::ReadyToDivide(std::vector<double> cellCycleInfluences)
{
	mpSimulationTime = SimulationTime::Instance();
	//assert(cellCycleInfluences.size()==0);
    bool ready;
        
    double timeSinceBirth = GetAge();
    
    switch (mCellType)
    {
        case STEM:
            ready = (timeSinceBirth >= mpCancerParams->GetStemCellCycleTime());
            break;
        case TRANSIT:
            ready = (timeSinceBirth >= RandomNumberGenerator::Instance()->NormalRandomDeviate(mpCancerParams->GetTransitCellCycleTime(), 1.0));
            break;
        default:
            ready = false;
            break;
    }
    
    return ready;
}
