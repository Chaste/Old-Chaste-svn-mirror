#include "TysonNovakCellCycleModel.hpp"

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

void TysonNovakCellCycleModel::ResetForDivision()
{	
    AbstractOdeBasedCellCycleModel::ResetForDivision();
    
    assert(mpOdeSystem!=NULL);
    
    /**
     * This model needs the protein concentrations and phase resetting to G0/G1.
     * 
     * In theory, the solution to the Tyson-Novak equations should exhibit stable 
     * oscillations, and we only need to halve the mass of the cell each period. 
     * 
     * However, the backward Euler solver used to solve the equations 
     * currently returns a solution that diverges after long times (see #316), so 
     * we must reset the initial conditions each period.
     */ 

    /// \todo: Uncomment this line and comment the line after once #316 is fixed
    /// mpOdeSystem->rGetStateVariables()[5] = mpOdeSystem->rGetStateVariables()[5]/2.0;

    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());
}

void TysonNovakCellCycleModel::InitialiseDaughterCell()
{    
    if (mpCell->GetCellType() == STEM) 
    {    
        mpCell->SetCellType(TRANSIT); 
    }
}

AbstractCellCycleModel* TysonNovakCellCycleModel::CreateDaughterCellCycleModel()
{
    return new TysonNovakCellCycleModel(mpOdeSystem->rGetStateVariables(), mBirthTime, mGeneration);
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

