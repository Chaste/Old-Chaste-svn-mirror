#include "AbstractCellCycleModel.hpp"

void AbstractCellCycleModel::SetCellType(CryptCellType cellType)
{
    mCellType = cellType;
}

void AbstractCellCycleModel::SetBirthTime(double birthTime)
{
    mBirthTime = birthTime;
}

double AbstractCellCycleModel::GetAge()
{
	mpSimulationTime = SimulationTime::Instance();
    return mpSimulationTime->GetDimensionalisedTime() - mBirthTime;
}

