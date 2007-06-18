#include "AbstractCellCycleModel.hpp"

// declare identifier for the serializer
// note that this has to be in the cpp file not the hpp
// Should only be needed for derived classes I think.
//BOOST_CLASS_EXPORT(AbstractCellCycleModel)

AbstractCellCycleModel::~AbstractCellCycleModel()
{}

void AbstractCellCycleModel::SetCellType(CryptCellType cellType)
{
    mCellType = cellType;
}

CryptCellType AbstractCellCycleModel::GetCellType()
{
    return mCellType;
}

CryptCellType AbstractCellCycleModel::UpdateCellType()
{   // This doesn't do anything in most classes, overwritten in others.
    return mCellType;
}

double AbstractCellCycleModel::GetBirthTime()
{
    return mBirthTime;
}

double AbstractCellCycleModel::GetAge()
{
    return SimulationTime::Instance()->GetDimensionalisedTime() - mBirthTime;
}

