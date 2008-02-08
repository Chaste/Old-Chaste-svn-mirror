#include "WntCellCycleModel.hpp"
#include "CellMutationStates.hpp"
#include "Exception.hpp"
#include <iostream>
#include <cassert>
#include <cfloat>

/**
 * A private constructor for daughter cells called by the CreateDaughterCellCycleModel function
 * (which can be called by TissueCell::CommonCopy() and isn't necessarily being born.
 *
 * @param pParentOdeSystem  to copy the state of.
 * @param rMutationState the mutation state of the cell (used by ODEs)
 * @param birthTime the simulation time when the cell divided (birth time of parent cell)
 * @param lastTime last time the cell cycle model was evaluated
 * @param inSG2MPhase whether the cell is in S-G2-M (not evaluating ODEs and just waiting)
 * @param readyToDivide 
 * @param divideTime If in the future this is the time at which the cell is going to divide
 */
WntCellCycleModel::WntCellCycleModel(AbstractOdeSystem* pParentOdeSystem,//const std::vector<double>& rParentProteinConcentrations,
                                     const CellMutationState& rMutationState, 
                                     double birthTime, double lastTime,
                                     bool inSG2MPhase, bool readyToDivide, double divideTime, unsigned generation)
   : AbstractWntOdeBasedCellCycleModel(lastTime) 
{
    if (pParentOdeSystem !=NULL)
    {
        std::vector<double> parent_protein_concs = pParentOdeSystem->rGetStateVariables();
        mpOdeSystem = new WntCellCycleOdeSystem(parent_protein_concs[8], rMutationState);// wnt pathway is reset in a couple of lines.
        // Set the model to be the same as the parent cell.
        mpOdeSystem->rGetStateVariables() = parent_protein_concs;
    }
    else
    {
        mpOdeSystem = NULL;
    }

    if (SimulationTime::Instance()->IsStartTimeSetUp()==false)
    {
        #define COVERAGE_IGNORE
        EXCEPTION("WntCellCycleModel is being created but SimulationTime has not been set up");
        #undef COVERAGE_IGNORE
    }
    mBirthTime = birthTime;
    mFinishedRunningOdes = inSG2MPhase;
    mReadyToDivide = readyToDivide;
    mDivideTime = divideTime;
    mGeneration = generation;
}

/**
 * A private constructor for archiving
 *
 * @param parentProteinConcentrations a std::vector of doubles of the protein concentrations (see WntCellCycleOdeSystem)
 * @param rMutationState the mutation state of the cell (used by ODEs)
 */
WntCellCycleModel::WntCellCycleModel(const std::vector<double>& rParentProteinConcentrations,
                                     const CellMutationState& rMutationState)
{
    mpOdeSystem = new WntCellCycleOdeSystem(rParentProteinConcentrations[8], rMutationState);// wnt pathway is reset in a couple of lines.
    // Set the model to be the same as the parent cell.
    mpOdeSystem->rGetStateVariables() = rParentProteinConcentrations;
}

/**
 * Returns a new WntCellCycleModel created with the correct initial conditions.
 *
 * Should be called just after the parent cell cycle model has been .Reset().
 */
AbstractCellCycleModel* WntCellCycleModel::CreateDaughterCellCycleModel()
{
    assert(mpCell!=NULL);
    // calls a cheeky version of the constructor which makes the new cell 
    // cycle model the same as the old one - not a dividing copy at this time.
    // unless the parent cell has just divided.
    return new WntCellCycleModel(mpOdeSystem, 
                                 mpCell->GetMutationState(), mBirthTime, mLastTime, 
                                 mFinishedRunningOdes, mReadyToDivide, mDivideTime, mGeneration);
}

/**
 * This carries out the work for ::UpdateCellType();
 * But does not check the current time so it can be used by the initialise method.
 */
void WntCellCycleModel::ChangeCellTypeDueToCurrentBetaCateninLevel()
{
    assert(mpOdeSystem!=NULL);
    assert(mpCell!=NULL);
    double beta_catenin_level = mpOdeSystem->rGetStateVariables()[6] + mpOdeSystem->rGetStateVariables()[7];
    //std::cout << "beta-catenin level = " << beta_catenin_level << "\n" << std::flush;        
    CellType cell_type=TRANSIT;
                
    // For mitogenic stimulus of 6x10^-4 in Wnt equations
    if (beta_catenin_level < 0.4127)
    {
        cell_type = DIFFERENTIATED;
    }
    
    mpCell->SetCellType(cell_type);
}

void WntCellCycleModel::Initialise()
{
    assert(mpOdeSystem==NULL);
    assert(mpCell!=NULL);
    mpOdeSystem = new WntCellCycleOdeSystem(WntConcentration::Instance()->GetWntLevel(mpCell), mpCell->GetMutationState());
    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());
    ChangeCellTypeDueToCurrentBetaCateninLevel();   
}    

bool WntCellCycleModel::SolveOdeToTime(double currentTime)
{
    // WE ARE IN G0 or G1 PHASE - running cell cycle ODEs
    double dt = 0.0001; // Needs to be this precise to stop crazy errors whilst we are still using rk4.

    // feed this time step's Wnt stimulus into the solver as a constant over this timestep.
    mpOdeSystem->rGetStateVariables()[8] = WntConcentration::Instance()->GetWntLevel(mpCell);
    // Use the cell's current mutation status as another input
    static_cast<WntCellCycleOdeSystem*>(mpOdeSystem)->SetMutationState(mpCell->GetMutationState());

    msSolver.SolveAndUpdateStateVariable(mpOdeSystem, mLastTime, currentTime, dt);

    mLastTime = currentTime;// normally done in Abstract class, but no harm in doing it here to prevent following line throwing an error.
    UpdateCellType();
    return msSolver.StoppingEventOccured();
}
    

