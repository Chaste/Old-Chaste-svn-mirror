#include "TysonNovakCellCycleModel.hpp"
#include "Exception.hpp"
#include <iostream>

TysonNovakCellCycleModel::TysonNovakCellCycleModel()
{
    mpSolver = new BackwardEulerIvpOdeSolver(6);
    mpSimulationTime = SimulationTime::Instance();
    if(mpSimulationTime->IsSimulationTimeSetUp()==false)
	{
        delete mpSolver;
		EXCEPTION("TysonNovakCellCycleModel is being created but SimulationTime has not been set up");
	}
    mLastTime = mpSimulationTime->GetDimensionalisedTime();
    mBirthTime = mLastTime;
    mDivideTime = mBirthTime;
    mOdeSystem.SetStateVariables(mOdeSystem.GetInitialConditions());
    mProteinConcentrations = mOdeSystem.GetInitialConditions();
    mReadyToDivide=false;
}

/**
 * A private constructor for daughter cells called only by the CreateCellCycleModel function
 * 
 * @param parentProteinConcentrations a std::vector of doubles of the protein concentrations (see WntCellCycleOdeSystem)
 * @param birthTime the SimulationTime when the cell divided (birth time of parent cell)
 */
TysonNovakCellCycleModel::TysonNovakCellCycleModel(std::vector<double> parentProteinConcentrations, double divideTime)
{
    mpSolver = new BackwardEulerIvpOdeSolver(6);
    mpSimulationTime = SimulationTime::Instance();
    if(mpSimulationTime->IsSimulationTimeSetUp()==false)
    {
        delete mpSolver;
        EXCEPTION("TysonNovakCellCycleModel is being created but SimulationTime has not been set up");
    }
    mLastTime = divideTime;
    mBirthTime = mLastTime;
    mDivideTime = mBirthTime;
    mProteinConcentrations = parentProteinConcentrations;
    mReadyToDivide=false;
}

TysonNovakCellCycleModel::~TysonNovakCellCycleModel()
{
    delete mpSolver;
}

void TysonNovakCellCycleModel::ResetModel()
{	// This model cycles itself and nothing needs to be reset.
	mBirthTime = mDivideTime;
    // std::cout << "divide time = " << mDivideTime << std::endl;
    // Halve the mass of the cell
    mProteinConcentrations[5] = mProteinConcentrations[5]/2.0;
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
    mDivideTime = birthTime;
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
        	double dt = 0.1/60.0;
            
    		OdeSolution solution = mpSolver->Solve(&mOdeSystem, mProteinConcentrations, mLastTime, current_time, dt, dt);
    	    
            unsigned timeRows = solution.GetNumberOfTimeSteps();
    	 	
    	 	for (unsigned i=0 ; i<6 ; i++)
    	 	{
    	 		mProteinConcentrations[i] = solution.rGetSolutions()[timeRows][i];
    	 		if (mProteinConcentrations[i]<0)
    	 		{
                    #define COVERAGE_IGNORE
    	 			std::cout << "Protein["<< i <<"] = "<< mProteinConcentrations[i] << "\n";
    	 			EXCEPTION("A protein concentration has gone negative\nCHASTE predicts that the TysonNovakCellCycleModel numerical method is probably unstable.");
                    #undef COVERAGE_IGNORE
    	 		}
    	 	}
    	 	
    	 	mLastTime = current_time;
    	    mReadyToDivide = mpSolver->StoppingEventOccured();
            if(mReadyToDivide)
            {
                unsigned end = solution.rGetSolutions().size() - 1;
                mDivideTime = solution.rGetTimes()[end];
            }
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
    return new TysonNovakCellCycleModel(mProteinConcentrations,mDivideTime);
}

