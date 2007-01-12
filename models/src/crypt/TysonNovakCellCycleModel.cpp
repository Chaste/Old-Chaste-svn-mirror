#include "TysonNovakCellCycleModel.hpp"
#include <iostream>

TysonNovakCellCycleModel::TysonNovakCellCycleModel()
{
    mpSimulationTime = SimulationTime::Instance();
    mLastTime = mpSimulationTime->GetDimensionalisedTime();
    
    mOdeSystem.SetStateVariables(mOdeSystem.GetInitialConditions());
    mProteinConcentrations = mOdeSystem.GetInitialConditions();
}
    
/// NOTE: the simulationTime parameter is NOT used!!!!!!!!!!!!!!
bool TysonNovakCellCycleModel::ReadyToDivide(double simulationTime)
{
    double current_time = mpSimulationTime->GetDimensionalisedTime();
    if(current_time<=mLastTime)
    {
    	EXCEPTION("TysonNovakCellCycleModel evaluated up to this time already");
    }
    
    double meshSize = 0.01;
	OdeSolution solution = mSolver.Solve(&mOdeSystem, mProteinConcentrations, mLastTime, current_time, meshSize, meshSize);

 	unsigned timeRows = solution.GetNumberOfTimeSteps();
 	
 	for (unsigned i=0 ; i<8 ; i++)
 	{
 		mProteinConcentrations[i] = solution.rGetSolutions()[timeRows][i];
 		if (mProteinConcentrations[i]<0)
 		{
 			std::cout << "Protein["<< i <<"] = "<< mProteinConcentrations[i] << "\n";
 			EXCEPTION("A protein concentration has gone negative\nCHASTE predicts that the TysonNovakCellCycleModel numerical method is probably unstable.");
 		}
 	}
 	
 	mLastTime = current_time;
    return mSolver.StoppingEventOccured();
}

std::vector<double> TysonNovakCellCycleModel::GetProteinConcentrations()
{
	return mProteinConcentrations;	
}
    
    
AbstractCellCycleModel* TysonNovakCellCycleModel::CreateCellCycleModel()
{
    return new TysonNovakCellCycleModel();
}

