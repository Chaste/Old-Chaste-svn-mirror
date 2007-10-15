#include "AbstractCellCycleModel.hpp"

AbstractCellCycleModel::~AbstractCellCycleModel()
{
    // Don't delete the cell - the cell deletes the cell cycle model 
    // when it is destroyed (not the following way round):
    // delete mpCell;
}

void AbstractCellCycleModel::SetCell(TissueCell* pCell)
{
    mpCell = pCell;
}

TissueCell* AbstractCellCycleModel::GetCell()
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

