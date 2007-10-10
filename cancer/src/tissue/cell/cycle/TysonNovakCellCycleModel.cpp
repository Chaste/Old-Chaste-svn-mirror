#include "TysonNovakCellCycleModel.hpp"
#include "Exception.hpp"
#include <iostream>

BackwardEulerIvpOdeSolver TysonNovakCellCycleModel::msSolver(6);

TysonNovakCellCycleModel::TysonNovakCellCycleModel()
{
    mpOdeSystem = new TysonNovak2001OdeSystem;
    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());
}

/**
 * A private constructor for daughter cells called only by the CreateCellCycleModel function
 *
 * @param parentProteinConcentrations a std::vector of doubles of the protein concentrations
 * @param birthTime the SimulationTime when the cell divided (birth time of parent cell)
 */
TysonNovakCellCycleModel::TysonNovakCellCycleModel(std::vector<double> parentProteinConcentrations, double divideTime)
 : AbstractOdeBasedCellCycleModel(divideTime)
{
    mpOdeSystem = new TysonNovak2001OdeSystem;
    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());
}

void TysonNovakCellCycleModel::ResetModel()
{	
    assert(mpOdeSystem!=NULL);
    // This model needs the protein concentrations and phase resetting to G0/G1.
    assert(mReadyToDivide);
    // This model should cycle itself and nothing needs to be reset.
    // but at the moment we are resetting to initial conditions because it
    // breaks after a while and will not converge.
    mBirthTime = mDivideTime;
    mLastTime = mDivideTime;
    
    //\TODO:Figure out why this goes unstable after a while...
    // Halve the mass of the cell
    if (false)
    {
        mpOdeSystem->rGetStateVariables()[5] = mpOdeSystem->rGetStateVariables()[5]/2.0;
    }
    else
    {
        mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());
    }

    mReadyToDivide = false;
}

bool TysonNovakCellCycleModel::ReadyToDivide()
{
    double current_time = SimulationTime::Instance()->GetDimensionalisedTime();
    
    if (!mReadyToDivide)
    {
        if (current_time>mLastTime)
        {
            double dt = 0.1/60.0;
            
            msSolver.SolveAndUpdateStateVariable(mpOdeSystem,mLastTime,current_time,dt);
            
            for (unsigned i=0 ; i<6 ; i++)
            {
                if (mpOdeSystem->rGetStateVariables()[i]<0)
                {
#define COVERAGE_IGNORE
                    std::cout << "Protein["<< i <<"] = "<< mpOdeSystem->rGetStateVariables()[i] << "\n";
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

AbstractCellCycleModel* TysonNovakCellCycleModel::CreateCellCycleModel()
{
    return new TysonNovakCellCycleModel(mpOdeSystem->rGetStateVariables(), mDivideTime);
}

