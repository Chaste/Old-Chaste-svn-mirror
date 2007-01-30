#include "TysonNovakCellCycleModel.hpp"
#include "Exception.hpp"
#include <iostream>

TysonNovakCellCycleModel::TysonNovakCellCycleModel()
{
    mpSimulationTime = SimulationTime::Instance();
    if(mpSimulationTime->IsSimulationTimeSetUp()==false)
	{
		EXCEPTION("TysonNovakCellCycleModel is being created but SimulationTime has not been set up");
	}
    mLastTime = mpSimulationTime->GetDimensionalisedTime();
    mBirthTime = mLastTime;
    mOdeSystem.SetStateVariables(mOdeSystem.GetInitialConditions());
    mProteinConcentrations = mOdeSystem.GetInitialConditions();
    mReadyToDivide=false;
}

void TysonNovakCellCycleModel::ResetModel()
{	// This model cycles itself and nothing needs to be reset.
	mBirthTime = mLastTime;
	mReadyToDivide=false;	
}

/**
 * Returns a new Tyson NovakCellCycleModel created with the correct initial conditions.
 * 
 * Should only be used in tests
 * 
 * @param birthTime the simulation time when the cell was born
 */  
void TysonNovakCellCycleModel::SetBirthTime(double birthTime)
{
	mLastTime = birthTime;
    mBirthTime = birthTime;
}
    
bool TysonNovakCellCycleModel::ReadyToDivide(std::vector<double> cellCycleInfluences)
{
	mpSimulationTime = SimulationTime::Instance();
	//assert(cellCycleInfluences.size()==0);
    double current_time = mpSimulationTime->GetDimensionalisedTime();
    
    if(!mReadyToDivide)
    {   
        if(current_time>mLastTime)
        {
        	double mesh_size = 0.1;
    		OdeSolution solution = mSolver.Solve(&mOdeSystem, mProteinConcentrations, mLastTime, current_time, mesh_size, mesh_size);
    	
            unsigned timeRows = solution.GetNumberOfTimeSteps();
    	 	
    	 	for (unsigned i=0 ; i<6 ; i++)
    	 	{
    	 		mProteinConcentrations[i] = solution.rGetSolutions()[timeRows][i];
    	 		if (mProteinConcentrations[i]<0)
    	 		{
    	 			std::cout << "Protein["<< i <<"] = "<< mProteinConcentrations[i] << "\n";
    	 			EXCEPTION("A protein concentration has gone negative\nCHASTE predicts that the TysonNovakCellCycleModel numerical method is probably unstable.");
    	 		}
    	 	}
    	 	
    	 	mLastTime = current_time;
    	    mReadyToDivide = mSolver.StoppingEventOccured();
        }
    }
    
    return mReadyToDivide;
}

std::vector<double> TysonNovakCellCycleModel::GetProteinConcentrations()
{
	return mProteinConcentrations;	
}
    
    
AbstractCellCycleModel* TysonNovakCellCycleModel::CreateCellCycleModel()
{
    return new TysonNovakCellCycleModel();
}

