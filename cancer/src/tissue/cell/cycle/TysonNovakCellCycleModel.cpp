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
 * A private constructor for daughter cells called only by the CreateDaughterCellCycleModel function
 *
 * @param parentProteinConcentrations a std::vector of doubles of the protein concentrations
 * @param birthTime the SimulationTime when the cell divided (birth time of parent cell)
 */
TysonNovakCellCycleModel::TysonNovakCellCycleModel(std::vector<double> parentProteinConcentrations, 
                                                   double divideTime, unsigned generation)
 : AbstractOdeBasedCellCycleModel(divideTime)
{
    mpOdeSystem = new TysonNovak2001OdeSystem;
    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());
    mGeneration = generation;
}

void TysonNovakCellCycleModel::ResetModel()
{	
    AbstractOdeBasedCellCycleModel::ResetModel();
    
    assert(mpOdeSystem!=NULL);
    // This model needs the protein concentrations and phase resetting to G0/G1.
    // This model should cycle itself and nothing needs to be reset.
    // but at the moment we are resetting to initial conditions because it
    // breaks after a while and will not converge.
    
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
    
}


AbstractCellCycleModel* TysonNovakCellCycleModel::CreateDaughterCellCycleModel()
{
    return new TysonNovakCellCycleModel(mpOdeSystem->rGetStateVariables(), mDivideTime, mGeneration);
}

bool TysonNovakCellCycleModel::SolveOdeToTime(double currentTime)
{
    double dt = 0.1/60.0;
    
    msSolver.SolveAndUpdateStateVariable(mpOdeSystem,mLastTime,currentTime,dt);
    
    return msSolver.StoppingEventOccured();
}

double TysonNovakCellCycleModel::GetOdeStopTime()
{
    assert(msSolver.StoppingEventOccured());
    return msSolver.GetStoppingTime();
}

/**
 * Tyson & Novak pretends it is running ODEs in just G1, 
 * but they really represent the whole cell cycle so 
 * we set the other phases to zero.
 */
double TysonNovakCellCycleModel::GetSDuration()
{
    return 0.0;
}

double TysonNovakCellCycleModel::GetG2Duration()
{
    return 0.0;
}

double TysonNovakCellCycleModel::GetMDuration()
{
    return 0.0;
}

