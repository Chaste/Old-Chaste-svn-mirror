#include "AbstractCellCycleModel.hpp"

// declare identifier for the serializer
// note that this has to be in the cpp file not the hpp
// Should only be needed for derived classes I think.
//BOOST_CLASS_EXPORT(AbstractCellCycleModel)

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

CryptCellType AbstractCellCycleModel::UpdateCellType()
{   // This doesn't do anything in most classes, overwritten in others.
    return mpCell->GetCellType();
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

