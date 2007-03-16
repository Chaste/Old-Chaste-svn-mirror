#include "AbstractCellCycleModel.hpp"

// declare identifier for the serializer
// note that this has to be in the cpp file not the hpp
// Should only be needed for derived classes I think.
//BOOST_CLASS_EXPORT(AbstractCellCycleModel)

AbstractCellCycleModel::~AbstractCellCycleModel()
{
}

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
    return SimulationTime::Instance()->GetDimensionalisedTime() - mBirthTime;
}

