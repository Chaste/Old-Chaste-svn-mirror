#include "WntCellCycleModel.hpp"
#include "Exception.hpp"
#include <iostream>
#include <cassert>

/**
 * Model constructor - an initial wnt stimulus must be provided to set up the
 * wnt pathway in an equilibrium state.
 *
 * Note that since we provide our own constructor, the compiler will *not*
 * generate a default one for us.
 * 
 * @param InitialWntStimulus a value between 0 and 1.
 * @param mutationStatus an unsigned taking the values 0 (healthy), 1 (APC+/-), 2 (BetaCat Delta45), 3 (APC-/-).
 *
 * \todo consider using an enum for the mutation state.
 */
WntCellCycleModel::WntCellCycleModel(double InitialWntStimulus, unsigned mutationStatus)
    : mOdeSystem(InitialWntStimulus, mutationStatus),
      mProteinConcentrations(mOdeSystem.GetInitialConditions())
{
    mpSimulationTime = SimulationTime::Instance();
    if(mpSimulationTime->IsStartTimeSetUp()==false)
    {
        EXCEPTION("WntCellCycleModel is being created but SimulationTime has not been set up");
    }
    mBirthTime = mpSimulationTime->GetDimensionalisedTime();
    mLastTime = mBirthTime;
    mInSG2MPhase = false;
    mReadyToDivide = false;
    mpCancerParams = CancerParameters::Instance();
}

/**
 * A private constructor for daughter cells called only by the CreateCellCycleModel function
 * 
 * @param parentProteinConcentrations a std::vector of doubles of the protein concentrations (see WntCellCycleOdeSystem)
 * @param birthTime the simulation time when the cell divided (birth time of parent cell)
 */
WntCellCycleModel::WntCellCycleModel(const std::vector<double>& rParentProteinConcentrations, double birthTime)
    : mOdeSystem(rParentProteinConcentrations[8], (unsigned)rParentProteinConcentrations[9]),
      mProteinConcentrations(mOdeSystem.GetInitialConditions())
{
    // Protein concentrations are initialised such that the cell cycle part of
    // the model is at the start of G1 phase.
	// Set the mutation state of the daughter to be the same as the parent cell.
	mProteinConcentrations[9] = rParentProteinConcentrations[9];
	// Set the Wnt pathway parts of the model to be the same as the parent cell.
	mProteinConcentrations[8] = rParentProteinConcentrations[8];
	mProteinConcentrations[7] = rParentProteinConcentrations[7];
	mProteinConcentrations[6] = rParentProteinConcentrations[6];
	mProteinConcentrations[5] = rParentProteinConcentrations[5];
    
    mpSimulationTime = SimulationTime::Instance();
    if(mpSimulationTime->IsStartTimeSetUp()==false)
	{
		EXCEPTION("WntCellCycleModel is being created but SimulationTime has not been set up");
	}
    mBirthTime = birthTime;
    mLastTime = birthTime;
    
    mInSG2MPhase = false;
    mReadyToDivide = false;
    
    mpCancerParams = CancerParameters::Instance();
}

WntCellCycleModel::~WntCellCycleModel()
{
}

/**
 * Resets the Wnt Model to the start of the cell cycle (this model does not cycle naturally)
 * Cells are given a new birth time and cell cycle proteins are reset.
 * Note that the wnt pathway proteins maintain their current values.
 * 
 * Should only be called by the MeinekeCryptCell Divide() method.
 * 
 */
void WntCellCycleModel::ResetModel()
{	// This model needs the protein concentrations and phase resetting to G0/G1.
	assert(mReadyToDivide);
	mLastTime = mDivideTime;
    mBirthTime = mDivideTime;
    // Keep the Wnt pathway in the same state but reset the cell cycle part
    // Cell cycle is proteins 0 to 4 (first 5 ODEs)
    for (unsigned i = 0 ; i<5 ; i++)
    {
	    mProteinConcentrations[i] = mOdeSystem.GetInitialConditions()[i];
    }
    mInSG2MPhase = false;
    mReadyToDivide = false;
}
    
/**
 * Returns a bool of whether or not the cell is ready to divide,
 * 
 * @param cellCycleInfluences the wnt stimulus -- must be the sole member of a standard vector of doubles and take a value between 0 and 1.
 * 
 */
bool WntCellCycleModel::ReadyToDivide(std::vector<double> cellCycleInfluences)
{
	//std::cout << "Looking up a cell cycle model" << std::endl;
	mpSimulationTime = SimulationTime::Instance();
	assert(cellCycleInfluences.size()==2);
	
 // Use the WntStimulus provided as an input
	mProteinConcentrations[8] = cellCycleInfluences[0];
 // Use the cell's current mutation status as another input
	mProteinConcentrations[9] = cellCycleInfluences[1];

	double current_time = mpSimulationTime->GetDimensionalisedTime();
	//std::cout << "Last time = " << mLastTime << ", Current Time = " << current_time << "\n" << std::endl;
	if(current_time>mLastTime)
	{
		if(!mInSG2MPhase)
		{	// WE ARE IN G0 or G1 PHASE - running cell cycle ODEs
			// feed this time step's Wnt stimulus into the solver as a constant over this timestep.
			
			
		    double meshSize = 0.0001; // Needs to be this precise to stop crazy errors whilst we are still using rk4.
			
//            for (unsigned i=0 ; i<8 ; i++)
//            
//                std::cout << "Before Protein["<< i <<"] = "<< mProteinConcentrations[i] << "\n";
//            }
//            std::cout.flush();
            OdeSolution solution = mSolver.Solve(&mOdeSystem, mProteinConcentrations, mLastTime, current_time, meshSize, meshSize);
//            std::cout << "After Solve\n";
//            std::cout.flush();
		
		 	unsigned timeRows = solution.GetNumberOfTimeSteps();
		 	if ( mSolver.StoppingEventOccured() == false )
		 	{
		 		//Check ODE solver for consistency
		 		assert (solution.rGetTimes()[timeRows] == current_time);
		 	}
		 	
		 	//std::cout<<"Stopping event = "<<mSolver.StoppingEventOccured()<<"\n";
		 	//std::cout<<"\n"<<last_time<<" == "<<current_time<<"\n";
		 	//
//            std::cout << "last time = "<<mLastTime<<"current time = "<<current_time<<"Number time steps = "<<timeRows<<"\n";
//            std::cout.flush();
		 	
		 	for (unsigned i=0 ; i<10 ; i++)
		 	{
		 		mProteinConcentrations[i] = solution.rGetSolutions()[timeRows][i];
//		 		std::cout << "Protein["<< i <<"] = "<< mProteinConcentrations[i] << "\n";
//		 		std::cout.flush();
                if (mProteinConcentrations[i]<0)
		 		{
                    #define COVERAGE_IGNORE
		 			std::cout << "Protein["<< i <<"] = "<< mProteinConcentrations[i] << "\n";
		 			EXCEPTION("A protein concentration has gone negative\nCHASTE predicts that the WntCellCycleModel numerical method is probably unstable.");
                    #undef COVERAGE_IGNORE
		 		}
		 	}
//            for (unsigned i=0 ; i<8 ; i++)
//            {
//                std::cout << "After Protein["<< i <<"] = "<< mProteinConcentrations[i] << "\n";
//            }
//            std::cout.flush();
						
			//std::cout << "Beta-Catenin = " << mProteinConcentrations[6] << "\n";
			if(mSolver.StoppingEventOccured())
		 	{
		 		unsigned end = solution.rGetSolutions().size() - 1;
    			// Tests the simulation is ending at the right time...(going into S phase at 5.971 hours)
    			double time_entering_S_phase = solution.rGetTimes()[end];
		 		mDivideTime = time_entering_S_phase + mpCancerParams->GetSG2MDuration();
		 		//std::cout << " Divide time = " << mDivideTime << "\n" << std::endl;
		 		mInSG2MPhase = true;
		 	}
//            std::cout << "After Stop check\n";
//            std::cout.flush();
		}
		else
		{	// WE ARE IN S-G2-M Phases, ODE model finished, just increasing time until division...
//	        std::cout << "In branch\n";
//            std::cout.flush();
    		if(current_time >= mDivideTime)
			{
				mReadyToDivide = true;
			}
		}
		mLastTime = current_time;
//    std::cout << "mLastTime updated to "<< mLastTime<<"\n";
//    std::cout.flush();
	}
//    std::cout << "Returning\n";
//    std::cout.flush();
    return mReadyToDivide;
}
  
  
/**
 * Returns the protein concentrations at the current time (useful for tests)
 * 
 */
std::vector<double> WntCellCycleModel::GetProteinConcentrations()
{
	return mProteinConcentrations;	
}

/**
 * Returns a new WntCellCycleModel created with the correct initial conditions.
 * 
 * Should be called just after the parent cell cycle model has been .Reset().
 * 
 */    
AbstractCellCycleModel* WntCellCycleModel::CreateCellCycleModel()
{
	// calls a cheeky version of the constructor which makes the new cell cycle model
	// the same age as the old one - not a copy at this time.
	return new WntCellCycleModel(mProteinConcentrations,mBirthTime);
}

/**
 * Returns a new WntCellCycleModel created with the correct initial conditions.
 * 
 * Should only be used in tests
 * 
 * @param birthTime the simulation time when the cell was born
 */  
void WntCellCycleModel::SetBirthTime(double birthTime)
{
	mLastTime = birthTime;
    mBirthTime = birthTime;
}

/**
 * Sets the protein concentrations and time when the model was last evaluated - should only be called by tests
 * 
 * @param lastTime the SimulationTime at which the protein concentrations apply
 * @param proteinConcentrations a standard vector of doubles of protein concentrations
 * 
 */  
void WntCellCycleModel::SetProteinConcentrationsForTestsOnly(double lastTime, std::vector<double> proteinConcentrations)
{
	mLastTime = lastTime;
	assert(proteinConcentrations.size()==10);
	mProteinConcentrations = proteinConcentrations;
}

