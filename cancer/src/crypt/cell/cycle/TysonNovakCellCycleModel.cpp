#include "TysonNovakCellCycleModel.hpp"
#include "Exception.hpp"
#include <iostream>

BackwardEulerIvpOdeSolver TysonNovakCellCycleModel::msSolver(6);

TysonNovakCellCycleModel::TysonNovakCellCycleModel()
    : AbstractCellCycleModel(),
      mOdeSystem()
{
    SimulationTime* p_sim_time = SimulationTime::Instance();
    if (p_sim_time->IsStartTimeSetUp()==false)
    {
        EXCEPTION("TysonNovakCellCycleModel is being created but SimulationTime has not been set up");
    }
    mLastTime = p_sim_time->GetDimensionalisedTime();
    mBirthTime = mLastTime;
    mDivideTime = mBirthTime;
    mReadyToDivide = false;
    mOdeSystem.SetStateVariables(mOdeSystem.GetInitialConditions());
}

/**
 * A private constructor for daughter cells called only by the CreateCellCycleModel function
 *
 * @param parentProteinConcentrations a std::vector of doubles of the protein concentrations
 * @param birthTime the SimulationTime when the cell divided (birth time of parent cell)
 */
TysonNovakCellCycleModel::TysonNovakCellCycleModel(std::vector<double> parentProteinConcentrations, double divideTime)
    : AbstractCellCycleModel(),
      mOdeSystem()
{
    if (SimulationTime::Instance()->IsStartTimeSetUp()==false)
    {
        EXCEPTION("TysonNovakCellCycleModel is being created but SimulationTime has not been set up");
    }
    mLastTime = divideTime;
    mBirthTime = divideTime;
    mDivideTime = divideTime;
    mReadyToDivide = false;
    
    mOdeSystem.SetStateVariables(parentProteinConcentrations);
}

TysonNovakCellCycleModel::~TysonNovakCellCycleModel()
{
}

void TysonNovakCellCycleModel::ResetModel()
{	// This model should cycle itself and nothing needs to be reset.
    // but at the moment we are resetting to initial conditions because it
    // breaks after a while and will not converge.
    mBirthTime = mDivideTime;
    mLastTime = mDivideTime;
    
    //\TODO:Figure out why this goes unstable after a while...
    // Halve the mass of the cell
    if (false)
    {
        mOdeSystem.rGetStateVariables()[5] = mOdeSystem.rGetStateVariables()[5]/2.0;
    }
    else
    {
        mOdeSystem.SetStateVariables(mOdeSystem.GetInitialConditions());
    }

    mReadyToDivide=false;
}

/**
 * Returns a new TysonNovakCellCycleModel created with the correct initial conditions.
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
    //assert(cellCycleInfluences.size()==0);
    double current_time = SimulationTime::Instance()->GetDimensionalisedTime();
    
    if (!mReadyToDivide)
    {
        if (current_time>mLastTime)
        {
            double dt = 0.1/60.0;
            
            msSolver.SolveAndUpdateStateVariable(&mOdeSystem,mLastTime,current_time,dt);
            
            for (unsigned i=0 ; i<6 ; i++)
            {
                if (mOdeSystem.rGetStateVariables()[i]<0)
                {
#define COVERAGE_IGNORE
                    std::cout << "Protein["<< i <<"] = "<< mOdeSystem.rGetStateVariables()[i] << "\n";
                    EXCEPTION("A protein concentration has gone negative\nCHASTE predicts that the TysonNovakCellCycleModel numerical method is probably unstable.");
#undef COVERAGE_IGNORE
                }
            }
            
            mLastTime = current_time;
            mReadyToDivide = msSolver.StoppingEventOccured();
            if (mReadyToDivide)
            {
                mDivideTime = msSolver.GetStoppingTime();
            }
        }
    }
    
    return mReadyToDivide;
}

std::vector<double> TysonNovakCellCycleModel::GetProteinConcentrations()
{
    return mOdeSystem.rGetStateVariables();
}


AbstractCellCycleModel* TysonNovakCellCycleModel::CreateCellCycleModel()
{
    return new TysonNovakCellCycleModel(mOdeSystem.rGetStateVariables(), mDivideTime);
}

