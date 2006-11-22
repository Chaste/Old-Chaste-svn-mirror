#include "TysonNovakCellCycleModel.hpp"
#include <iostream>

TysonNovakCellCycleModel::TysonNovakCellCycleModel()
{
    mpSimulationTime = SimulationTime::Instance();
    mLastTime = mpSimulationTime->GetDimensionalisedTime();
    
    mOdeSystem.SetStateVariables(mOdeSystem.GetInitialConditions());    
}
    
/// NOTE: the simulationTime parameter is NOT used!!!!!!!!!!!!!!
bool TysonNovakCellCycleModel::ReadyToDivide(double simulationTime)
{
    double current_time = mpSimulationTime->GetDimensionalisedTime();

    std::vector<double> state_variables = mOdeSystem.rGetStateVariables();

    mSolver.Solve(&mOdeSystem, state_variables, mLastTime, current_time, 0.01, 0.01);
 
    return mSolver.StoppingEventOccured();
}
    
    
AbstractCellCycleModel* TysonNovakCellCycleModel::CreateCellCycleModel()
{
    return new TysonNovakCellCycleModel();
}

