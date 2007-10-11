#include "WntCellCycleModel.hpp"
#include "CellMutationStates.hpp"
#include "Exception.hpp"
#include "WntGradient.hpp"
#include <iostream>
#include <cassert>
#include <cfloat>

RungeKutta4IvpOdeSolver WntCellCycleModel::msSolver;

/**
 * Model constructor
 *
 * Note that since we provide our own constructor, the compiler will *not*
 * generate a default one for us. *
 */
WntCellCycleModel::WntCellCycleModel()
{
    mInSG2MPhase = false;
}

/**
 * A private constructor for daughter cells called by the CreateCellCycleModel function
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
                                     bool inSG2MPhase, bool readyToDivide, double divideTime)
   : AbstractOdeBasedCellCycleModel(lastTime) 
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
    SimulationTime* p_sim_time = SimulationTime::Instance();
    if (p_sim_time->IsStartTimeSetUp()==false)
    {
        #define COVERAGE_IGNORE
        EXCEPTION("WntCellCycleModel is being created but SimulationTime has not been set up");
        #undef COVERAGE_IGNORE
    }
    mBirthTime = birthTime;
    mInSG2MPhase = inSG2MPhase;
    mReadyToDivide = readyToDivide;
    mDivideTime = divideTime;
}

/**
 * A private constructor for archiving
 *
 * @param parentProteinConcentrations a std::vector of doubles of the protein concentrations (see WntCellCycleOdeSystem)
 * @param rMutationState the mutation state of the cell (used by ODEs)
 * @param birthTime the simulation time when the cell divided (birth time of parent cell)
 * @param lastTime last time the cell cycle model was evaluated
 * @param inSG2MPhase whether the cell is in S-G2-M (not evaluating ODEs and just waiting)
 * @param readyToDivide 
 * @param divideTime If in the future this is the time at which the cell is going to divide
 */
WntCellCycleModel::WntCellCycleModel(const std::vector<double>& rParentProteinConcentrations,
                                     const CellMutationState& rMutationState, 
                                     double birthTime, double lastTime,
                                     bool inSG2MPhase, bool readyToDivide, double divideTime)
    : AbstractOdeBasedCellCycleModel(lastTime) 
{
    mpOdeSystem = new WntCellCycleOdeSystem(rParentProteinConcentrations[8], rMutationState);// wnt pathway is reset in a couple of lines.
    // Set the model to be the same as the parent cell.
    mpOdeSystem->rGetStateVariables() = rParentProteinConcentrations;
    
    SimulationTime* p_sim_time = SimulationTime::Instance();
    if (p_sim_time->IsStartTimeSetUp()==false)
    {
        #define COVERAGE_IGNORE
        EXCEPTION("WntCellCycleModel is being created but SimulationTime has not been set up");
        #undef COVERAGE_IGNORE
    }
    mBirthTime = birthTime;
    mInSG2MPhase = inSG2MPhase;
    mReadyToDivide = readyToDivide;
    mDivideTime = divideTime;
}

/**
 * Resets the Wnt Model to the start of the cell cycle (this model does not cycle naturally)
 * Cells are given a new birth time and cell cycle proteins are reset.
 * Note that the wnt pathway proteins maintain their current values.
 *
 * Should only be called by the TissueCell::Divide() method.
 */
void WntCellCycleModel::ResetModel()
{	
    assert(mpOdeSystem!=NULL);
    
    // This model needs the protein concentrations and phase resetting to G0/G1.
    assert(mReadyToDivide);
    mLastTime = mDivideTime;
    mBirthTime = mDivideTime;
    // Keep the Wnt pathway in the same state but reset the cell cycle part
    // Cell cycle is proteins 0 to 4 (first 5 ODEs)
    std::vector<double> init_conds = mpOdeSystem->GetInitialConditions();
    for (unsigned i = 0 ; i<5 ; i++)
    {
        mpOdeSystem->rGetStateVariables()[i] = init_conds[i];
    }
    mInSG2MPhase = false;
    mReadyToDivide = false;
}

/**
 * Returns a bool of whether or not the cell is ready to divide,
 *
 * @param cellCycleInfluences the wnt stimulus -- must be the sole member 
 * of a standard vector of doubles and take a value between 0 and 1.
 * 
 * This function sets one of the ODE system variables to represent the 
 * Wnt level and also gives the ODE system the mutation state of this cell.
 */
bool WntCellCycleModel::ReadyToDivide()
{
    assert(mpOdeSystem!=NULL);
    assert(mpCell!=NULL);
    
    double current_time = SimulationTime::Instance()->GetDimensionalisedTime();
    
    if (current_time>mLastTime)
    {
        if (!mInSG2MPhase)
        {	// WE ARE IN G0 or G1 PHASE - running cell cycle ODEs
            double dt = 0.0001; // Needs to be this precise to stop crazy errors whilst we are still using rk4.

            // feed this time step's Wnt stimulus into the solver as a constant over this timestep.
            mpOdeSystem->rGetStateVariables()[8] = WntGradient::Instance()->GetWntLevel(mpCell);
            // Use the cell's current mutation status as another input
            static_cast<WntCellCycleOdeSystem*>(mpOdeSystem)->SetMutationState(mpCell->GetMutationState());
    
            msSolver.SolveAndUpdateStateVariable(mpOdeSystem, mLastTime, current_time, dt);

            for (unsigned i=0 ; i<mpOdeSystem->GetNumberOfStateVariables() ; i++)
            {
                if (mpOdeSystem->rGetStateVariables()[i]<0)
                {
                    #define COVERAGE_IGNORE
                    std::cout << "Protein["<< i <<"] = "<< mpOdeSystem->rGetStateVariables()[i] << "\n";
                    EXCEPTION("A protein concentration has gone negative\nCHASTE predicts that the WntCellCycleModel numerical method is probably unstable.");
                    #undef COVERAGE_IGNORE
                }
            }

            if (msSolver.StoppingEventOccured())
            {
                // Tests the simulation is ending at the right time...(going into S phase at 5.971 hours)
                double time_entering_S_phase = msSolver.GetStoppingTime();
                mDivideTime = time_entering_S_phase + GetWntSG2MDuration();

                mInSG2MPhase = true;
            }
            mLastTime = current_time;   // This is the last time the ODEs were evaluated.
            UpdateCellType();
        }
        else
        {	
            // WE ARE IN S-G2-M Phases, ODE model finished, just increasing time until division...
            if (current_time >= mDivideTime)
            {
                mReadyToDivide = true;
            }
        }
    }
    
    return mReadyToDivide;
}

/**
 * Returns a new WntCellCycleModel created with the correct initial conditions.
 *
 * Should be called just after the parent cell cycle model has been .Reset().
 */
AbstractCellCycleModel* WntCellCycleModel::CreateCellCycleModel()
{
    assert(mpCell!=NULL);
    // calls a cheeky version of the constructor which makes the new cell 
    // cycle model the same as the old one - not a dividing copy at this time.
    // unless the parent cell has just divided.
    return new WntCellCycleModel(mpOdeSystem, 
                                 mpCell->GetMutationState(), mBirthTime, mLastTime, 
                                 mInSG2MPhase, mReadyToDivide, mDivideTime);
}

/**
 * Updates the current cell type to reflect whether the
 * beta-catenin level has dropped low enough to make it stop dividing.
 * This should only be called when the cell cycle model has been 
 * evaluated to the current time, or it may give misleading results.
 */
void WntCellCycleModel::UpdateCellType()
{
    assert(mpOdeSystem!=NULL);
    assert(mpCell!=NULL);
    if (SimulationTime::Instance()->GetDimensionalisedTime() > mLastTime)
    {
        EXCEPTION("WntCellCycleModel::UpdateCellType() should only be called when the cell cycle model has been evaluated to the current time\n");   
    }
    ChangeCellTypeDueToCurrentBetaCateninLevel();
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

double WntCellCycleModel::GetWntSG2MDuration()
{   
    // overridden in subclass StochasticWntCellCycleModel
    return CancerParameters::Instance()->GetSG2MDuration();
}

void WntCellCycleModel::Initialise()
{
    assert(mpOdeSystem==NULL);
    assert(mpCell!=NULL);
    mpOdeSystem = new WntCellCycleOdeSystem(WntGradient::Instance()->GetWntLevel(mpCell), mpCell->GetMutationState());
    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());
    ChangeCellTypeDueToCurrentBetaCateninLevel();   
}    

