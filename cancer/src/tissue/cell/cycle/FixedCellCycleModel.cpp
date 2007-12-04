#include "FixedCellCycleModel.hpp"

AbstractCellCycleModel *FixedCellCycleModel::CreateCellCycleModel()
{
    return new FixedCellCycleModel(mG1Duration, mGeneration);
}

