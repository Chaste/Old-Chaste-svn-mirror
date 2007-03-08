#include "StochasticCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "Exception.hpp"
#include <cassert>

// For some reason(??) the archiver doesn't want the following lines for this class!

// declare identifier for the serializer
// note that this has to be in the cpp file not the hpp
// BOOST_CLASS_EXPORT_GUID(StochasticCellCycleModel, "StochasticCellCycleModel")

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
	mpCancerParams = CancerParameters::Instance();
	mBirthTime = mpSimulationTime->GetDimensionalisedTime();
    mpGen=pGen;
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
            ready = (timeSinceBirth >= mpGen->NormalRandomDeviate(mpCancerParams->GetTransitCellCycleTime(), 1.0));
            break;
        default:
            ready = false;
            break;
    }
    
    return ready;
}
