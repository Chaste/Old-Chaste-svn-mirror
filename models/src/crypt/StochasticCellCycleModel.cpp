#include "StochasticCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "Exception.hpp"


AbstractCellCycleModel *StochasticCellCycleModel::CreateCellCycleModel()
{
    return new StochasticCellCycleModel(mpGen);
}


StochasticCellCycleModel::StochasticCellCycleModel(RandomNumberGenerator *pGen)
{
	mpSimulationTime = SimulationTime::Instance();
	if(mpSimulationTime->IsSimulationTimeSetUp()==false)
	{
		EXCEPTION("StochasticCellCycleModel is being created but SimulationTime has not been set up");
	}
	mp_params = CancerParameters::Instance();
	mBirthTime = mpSimulationTime->GetDimensionalisedTime();
    mpGen=pGen;
}

void StochasticCellCycleModel::ResetModel()
{
	mpSimulationTime = SimulationTime::Instance();
	mBirthTime = mpSimulationTime->GetDimensionalisedTime();	
}

bool StochasticCellCycleModel::ReadyToDivide(std::vector<double> cellCycleInfluences)
{
	mpSimulationTime = SimulationTime::Instance();
	assert(cellCycleInfluences.size()==0);
    bool ready;
        
    double timeSinceBirth = GetAge();
    
    switch (mCellType)
    {
        case STEM:
            ready = (timeSinceBirth >= mp_params->GetStemCellCycleTime());
            break;
        case TRANSIT:
            ready = (timeSinceBirth >= mpGen->NormalRandomDeviate(mp_params->GetTransitCellCycleTime(), 1.0));
            break;
        default:
            ready = false;
            break;
    }
    
    return ready;
}
