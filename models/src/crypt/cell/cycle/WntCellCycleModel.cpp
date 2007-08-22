#include "WntCellCycleModel.hpp"
#include "WntGradient.hpp"
#include "CryptCellMutationStates.hpp"
#include "Exception.hpp"
#include <iostream>
#include <cassert>
#include <cfloat>

RungeKutta4IvpOdeSolver WntCellCycleModel::msSolver;

/**
 * Model constructor - an initial wnt stimulus must be provided to set up the
 * wnt pathway in an equilibrium state.
 *
 * Note that since we provide our own constructor, the compiler will *not*
 * generate a default one for us.
 *
 * @param InitialWntStimulus a value between 0 and 1.
 */
WntCellCycleModel::WntCellCycleModel(double InitialWntStimulus, WntGradient &rWntGradient)
        : AbstractCellCycleModel(),
          mpOdeSystem(NULL),
          mInitialWntStimulus(InitialWntStimulus),
          mrWntGradient(rWntGradient)          
{
    SimulationTime* p_sim_time = SimulationTime::Instance();
    if (p_sim_time->IsStartTimeSetUp()==false)
    {
        EXCEPTION("WntCellCycleModel is being created but SimulationTime has not been set up");
    }
    mBirthTime = p_sim_time->GetDimensionalisedTime();
    mLastTime = mBirthTime;
    mInSG2MPhase = false;
    mReadyToDivide = false;
    mDivideTime = DBL_MAX;
}

/**
 * A private constructor for daughter cells called only by the CreateCellCycleModel function
 *
 * @param parentProteinConcentrations a std::vector of doubles of the protein concentrations (see WntCellCycleOdeSystem)
 * @param birthTime the simulation time when the cell divided (birth time of parent cell)
 */
WntCellCycleModel::WntCellCycleModel(const std::vector<double>& rParentProteinConcentrations,
                                     const CryptCellMutationState& rMutationState, 
                                     double birthTime, double lastTime, 
                                     WntGradient& rWntGradient,
                                     bool inSG2MPhase, bool readyToDivide, double divideTime)
        : AbstractCellCycleModel(),
          mrWntGradient(rWntGradient)
{
    mpOdeSystem = new WntCellCycleOdeSystem(rParentProteinConcentrations[8], rMutationState);// wnt pathway is reset in a couple of lines.
    // Set the model to be the same as the parent cell.
    mpOdeSystem->rGetStateVariables() = rParentProteinConcentrations;
    
    SimulationTime* p_sim_time = SimulationTime::Instance();
    if (p_sim_time->IsStartTimeSetUp()==false)
    {
        EXCEPTION("WntCellCycleModel is being created but SimulationTime has not been set up");
    }
    mBirthTime = birthTime;
    mLastTime = lastTime;
    mInSG2MPhase = inSG2MPhase;
    mReadyToDivide = readyToDivide;
    mDivideTime = divideTime;
}

WntCellCycleModel::~WntCellCycleModel()
{
    delete mpOdeSystem;
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
bool WntCellCycleModel::ReadyToDivide(std::vector<double> cellCycleInfluences)
{
    assert(mpOdeSystem!=NULL);
    assert(mpCell!=NULL);
    assert(cellCycleInfluences.size()==1);
    
    // Use the WntStimulus provided as an input
    mpOdeSystem->rGetStateVariables()[8] = cellCycleInfluences[0];
    // Use the cell's current mutation status as another input
    mpOdeSystem->SetMutationState(mpCell->GetMutationState());
    
    double current_time = SimulationTime::Instance()->GetDimensionalisedTime();
    
    if (current_time>mLastTime)
    {
        if (!mInSG2MPhase)
        {	// WE ARE IN G0 or G1 PHASE - running cell cycle ODEs
            // feed this time step's Wnt stimulus into the solver as a constant over this timestep.
            double dt = 0.0001; // Needs to be this precise to stop crazy errors whilst we are still using rk4.
            
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
        }
        else
        {	
            // WE ARE IN S-G2-M Phases, ODE model finished, just increasing time until division...
            if (current_time >= mDivideTime)
            {
                mReadyToDivide = true;
            }
        }
        mLastTime = current_time;
    }
    return mReadyToDivide;
}


/**
 * Returns the protein concentrations at the current time (useful for tests)
 *
 * NB: Will copy the vector - you can't use this to modify the concentrations.
 */
std::vector<double>  WntCellCycleModel::GetProteinConcentrations() const
{
    assert(mpOdeSystem!=NULL);
    return mpOdeSystem->rGetStateVariables();
}

/**
 * Returns a new WntCellCycleModel created with the correct initial conditions.
 *
 * Should be called just after the parent cell cycle model has been .Reset().
 */
AbstractCellCycleModel* WntCellCycleModel::CreateCellCycleModel()
{
    assert(mpOdeSystem!=NULL);
    assert(mpCell!=NULL);
    // calls a cheeky version of the constructor which makes the new cell 
    // cycle model the same as the old one - not a dividing copy at this time.
    // unless the parent cell has just divided.
    return new WntCellCycleModel(mpOdeSystem->rGetStateVariables(), 
                                mpCell->GetMutationState(), mBirthTime, mLastTime, 
                                mrWntGradient, mInSG2MPhase, mReadyToDivide, mDivideTime);
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
    assert(mpOdeSystem!=NULL);
    assert(proteinConcentrations.size()==mpOdeSystem->rGetStateVariables().size());
    mLastTime = lastTime;
    mpOdeSystem->SetStateVariables(proteinConcentrations);
}


CryptCellType WntCellCycleModel::UpdateCellType()
{
    assert(mpOdeSystem!=NULL);
    double betaCateninLevel = mpOdeSystem->rGetStateVariables()[6] + mpOdeSystem->rGetStateVariables()[7];
    //std::cout << "beta-catenin level = " << betaCateninLevel << "\n" << std::flush;        
    CryptCellType cell_type=TRANSIT;
                
    // For mitogenic stimulus of 6x10^-4 in Wnt equations
    if (betaCateninLevel < 0.4127)
    {
        cell_type = DIFFERENTIATED;
    }
    return cell_type;
}

double WntCellCycleModel::GetWntSG2MDuration()
{   // overridden in subclass StochasticWntCellCycleModel
    return CancerParameters::Instance()->GetSG2MDuration();
}

void WntCellCycleModel::SetCell(MeinekeCryptCell* pCell)
{   
    mpCell = pCell;
    if(mpOdeSystem==NULL)   // this is a new cell (not a dividing original cell) and needs an ODE system.
    { 
        mpOdeSystem = new WntCellCycleOdeSystem(mInitialWntStimulus, pCell->GetMutationState());
        mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());
    }            
}
