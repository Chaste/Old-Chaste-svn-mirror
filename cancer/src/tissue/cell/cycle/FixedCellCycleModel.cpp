#include "FixedCellCycleModel.hpp"

AbstractCellCycleModel *FixedCellCycleModel::CreateDaughterCellCycleModel()
{
    return new FixedCellCycleModel(mG1Duration, mGeneration);
}

