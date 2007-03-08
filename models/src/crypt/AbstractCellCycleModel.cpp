#include "AbstractCellCycleModel.hpp"

// declare identifier for the serializer
// note that this has to be in the cpp file not the hpp
BOOST_CLASS_EXPORT_GUID(AbstractCellCycleModel, "AbstractCellCycleModel")


void AbstractCellCycleModel::SetCellType(CryptCellType cellType)
{
    mCellType = cellType;
}

CryptCellType AbstractCellCycleModel::GetCellType()
{
    return mCellType;
}

//void AbstractCellCycleModel::SetBirthTime(double birthTime)
//{
//    mBirthTime = birthTime;
//}

double AbstractCellCycleModel::GetBirthTime()
{
	return mBirthTime;
}

double AbstractCellCycleModel::GetAge()
{
	mpSimulationTime = SimulationTime::Instance();
    return mpSimulationTime->GetDimensionalisedTime() - mBirthTime;
}

