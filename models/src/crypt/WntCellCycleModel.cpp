#include "WntCellCycleModel.hpp"
#include "Exception.hpp"
#include <iostream>
#include <cassert>

WntCellCycleModel::WntCellCycleModel()
{
	EXCEPTION("A Wnt cell cycle model must be given a steady state of Wnt (double)\n to set steady state of wnt pathway at start of model");	
}

WntCellCycleModel::WntCellCycleModel(double InitialWntStimulus)
{
	WntCellCycleOdeSystem mOdeSystem(InitialWntStimulus);
    mpSimulationTime = SimulationTime::Instance();
    if(mpSimulationTime->IsSimulationTimeSetUp()==false)
	{
		EXCEPTION("WntCellCycleModel is being created but SimulationTime has not been set up");
	}
    mLastTime = mpSimulationTime->GetDimensionalisedTime();
    mBirthTime = mLastTime;
    mProteinConcentrations = mOdeSystem.GetInitialConditions();
    mInSG2MPhase = false;
    mpCancerParams = CancerParameters::Instance();
}

void WntCellCycleModel::ResetModel()
{	// This model needs the protein concentrations and phase resetting to G0/G1.
	mpSimulationTime = SimulationTime::Instance();
	mLastTime = mpSimulationTime->GetDimensionalisedTime();
    mBirthTime = mLastTime;
    // Keep the Wnt pathway in the same state but reset the cell cycle part
    // Cell cycle is proteins 0 to 4 (first 5 ODEs)
    for (unsigned i = 0 ; i<5 ; i++)
    {
	    mProteinConcentrations[i] = mOdeSystem.GetInitialConditions()[i];
    }
    mInSG2MPhase = false;
}
    
bool WntCellCycleModel::ReadyToDivide(std::vector<double> cellCycleInfluences)
{
	mpSimulationTime = SimulationTime::Instance();
	assert(cellCycleInfluences.size()==1);
	double wntStimulus = cellCycleInfluences[0];
	bool divideNow = false;
	double current_time = mpSimulationTime->GetDimensionalisedTime();
	
	if(!mInSG2MPhase)
	{	// WE ARE IN G0 or G1 PHASE - running cell cycle ODEs
		// feed this time step's Wnt stimulus into the solver as a constant over this timestep.
		mProteinConcentrations[7] = wntStimulus;
		
		if(current_time<=mLastTime)
	    {
	    	EXCEPTION("WntCellCycleModel evaluated up to this time already");
	    }
	    
		double meshSize = 0.0001; // Needs to be this precise to stop crazy errors whilst we are still using rk4.
		OdeSolution solution = mSolver.Solve(&mOdeSystem, mProteinConcentrations, mLastTime, current_time, meshSize, meshSize);
	
	 	unsigned timeRows = solution.GetNumberOfTimeSteps();
	 	
	 	for (unsigned i=0 ; i<8 ; i++)
	 	{
	 		mProteinConcentrations[i] = solution.rGetSolutions()[timeRows][i];
	 		if (mProteinConcentrations[i]<0)
	 		{
	 			std::cout << "Protein["<< i <<"] = "<< mProteinConcentrations[i] << "\n";
	 			EXCEPTION("A protein concentration has gone negative\nCHASTE predicts that the WntCellCycleModel numerical method is probably unstable.");
	 		}
	 	}
		mLastTime = current_time;
		
		//std::cout << "Beta-Catenin = " << mProteinConcentrations[6] << "\n";
		if(mSolver.StoppingEventOccured())
	 	{
	 		mDivideTime = current_time + mpCancerParams->GetSG2MDuration();
	 		mInSG2MPhase = true;
	 	}
	}
	else
	{	// WE ARE IN S-G2-M Phases, ODE model finished just increasing time until division...
		if(current_time >= mDivideTime)
		{
			divideNow = true;
		}
	}
 	
    return divideNow;
}
    
std::vector<double> WntCellCycleModel::GetProteinConcentrations()
{
	return mProteinConcentrations;	
}
    
AbstractCellCycleModel* WntCellCycleModel::CreateCellCycleModel()
{
    return new WntCellCycleModel();
}

