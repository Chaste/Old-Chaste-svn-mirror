#include "AbstractSimpleMeinekeCellCycleModel.hpp"

void AbstractSimpleMeinekeCellCycleModel::ResetModel()
{ 
//    
//    if( mpCell->GetCellType() == TRANSIT)
//    {
//        mGeneration++;
//    }
    mBirthTime = SimulationTime::Instance()->GetDimensionalisedTime();
    SetG1Duration();
}

