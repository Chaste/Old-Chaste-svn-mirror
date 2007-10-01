#include "AbstractCellCycleModel.hpp"

AbstractCellCycleModel::~AbstractCellCycleModel()
{
    //delete mpCell;
}

void AbstractCellCycleModel::SetCell(MeinekeCryptCell* pCell)
{
    mpCell=pCell;
}

MeinekeCryptCell* AbstractCellCycleModel::GetCell()
{
    assert(mpCell!=NULL);
    return mpCell;
}

void AbstractCellCycleModel::SetBirthTime(double birthTime)
{
    mBirthTime = birthTime;
}

double AbstractCellCycleModel::GetBirthTime() const
{
    return mBirthTime;
}

double AbstractCellCycleModel::GetAge()
{
    return SimulationTime::Instance()->GetDimensionalisedTime() - mBirthTime;
}

