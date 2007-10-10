#include "OxygenBasedCellCycleModel.hpp"
#include "Exception.hpp"
#include <iostream>
#include <cassert>
#include <cfloat>

RungeKutta4IvpOdeSolver OxygenBasedCellCycleModel::msSolver;

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
OxygenBasedCellCycleModel::OxygenBasedCellCycleModel(AbstractOdeSystem* pParentOdeSystem, 
                              const CellMutationState& rMutationState, double birthTime, 
                              double lastTime, bool readyToDivide, double divideTime)
    : AbstractOdeBasedCellCycleModel(lastTime)// these values overwritten below
{
    if (pParentOdeSystem !=NULL)
    {
        std::vector<double> parent_protein_concs = pParentOdeSystem->rGetStateVariables();
        mpOdeSystem = new Alarcon2004OxygenBasedCellCycleOdeSystem(parent_protein_concs[5], rMutationState);
        
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
        EXCEPTION("OxygenBasedCellCycleModel is being created but SimulationTime has not been set up");
        #undef COVERAGE_IGNORE
    }
    mBirthTime = birthTime;
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
OxygenBasedCellCycleModel::OxygenBasedCellCycleModel(const std::vector<double>& rParentProteinConcentrations, 
                              const CellMutationState& rMutationState, double birthTime, 
                              double lastTime, bool readyToDivide, double divideTime)
    : AbstractOdeBasedCellCycleModel(lastTime)// these values overwritten below
{
    mpOdeSystem = new Alarcon2004OxygenBasedCellCycleOdeSystem(rParentProteinConcentrations[5], rMutationState);
    
    // Set the model to be the same as the parent cell.
    mpOdeSystem->rGetStateVariables() = rParentProteinConcentrations;
    
    SimulationTime* p_sim_time = SimulationTime::Instance();
    if (p_sim_time->IsStartTimeSetUp()==false)
    {
        #define COVERAGE_IGNORE
        EXCEPTION("OxygenBasedCellCycleModel is being created but SimulationTime has not been set up");
        #undef COVERAGE_IGNORE
    }
    mBirthTime = birthTime;
    mReadyToDivide = readyToDivide;
    mDivideTime = divideTime;
}

/**
 * Resets the oxygen-based model to the start of the cell cycle (this model does not cycle naturally)
 * Cells are given a new birth time and cell cycle proteins are reset.
 * Note that the oxygen concentration maintains its current value.
 *
 * Should only be called by the TissueCell Divide() method.
 *
 */
void OxygenBasedCellCycleModel::ResetModel()
{	
    assert(mpOdeSystem!=NULL);
    
    // This model needs the protein concentrations and phase resetting to G0/G1.
    assert(mReadyToDivide);
    mLastTime = mDivideTime;
    
    mBirthTime = mDivideTime;
    // Keep the oxygen concentration the same but reset everything else
    std::vector<double> init_conds = mpOdeSystem->GetInitialConditions();
    for (unsigned i = 0 ; i<5 ; i++)
    {
        mpOdeSystem->rGetStateVariables()[i] = init_conds[i];
    }
    mReadyToDivide = false;
}

/**
 * Returns a bool of whether or not the cell is ready to divide.
 *
 * This method sets one of the ODE system variables to represent the 
 * oxygen concentration and gives the ODE system the mutation state of this cell.
 */
bool OxygenBasedCellCycleModel::ReadyToDivide()
{
    assert(mpOdeSystem!=NULL);
    assert(mpCell!=NULL);
    
    // Danger Will Robinson! DIM currently hard-coded to 2
    mpOdeSystem->rGetStateVariables()[5] = CellwiseData<2>::Instance()->GetValue(mpCell,0);
     
    // Use the cell's current mutation status as another input
    static_cast<Alarcon2004OxygenBasedCellCycleOdeSystem*>(mpOdeSystem)->SetMutationState(mpCell->GetMutationState());
    
    double current_time = SimulationTime::Instance()->GetDimensionalisedTime();
    
    if (current_time>mLastTime)
    {
        // feed this time step's oxygen concentration into the solver as a constant over this timestep.
        double dt = 0.001; // Needs to be this precise to stop crazy errors whilst we are still using rk4.
        msSolver.SolveAndUpdateStateVariable(mpOdeSystem, mLastTime, current_time, dt);

        // the exception for oxygen concentration < 0 is currently handled elsewhere            
        for (unsigned i=0 ; i<mpOdeSystem->GetNumberOfStateVariables()-1 ; i++)
        {
            if (mpOdeSystem->rGetStateVariables()[i]<0)
            {
                #define COVERAGE_IGNORE
                std::cout << "Protein["<< i <<"] = "<< mpOdeSystem->rGetStateVariables()[i] << "\n";
                EXCEPTION("A protein concentration has gone negative\nCHASTE predicts that the OxygenBasedCellCycleModel numerical method is probably unstable.");
                #undef COVERAGE_IGNORE
            }
        }

        if (msSolver.StoppingEventOccured())
        {
            mReadyToDivide = true;
            
            // Need to test that the simulation is ending at the right time
            double time_entering_S_phase = msSolver.GetStoppingTime();
            mDivideTime = time_entering_S_phase + CancerParameters::Instance()->GetSG2MDuration();

        }
        mLastTime = current_time;
    }
    
    return mReadyToDivide;
}


/**
 * Returns a new OxygenBasedCellCycleModel created with the correct initial conditions.
 *
 * Should be called just after the parent cell cycle model has been .Reset().
 */
AbstractCellCycleModel* OxygenBasedCellCycleModel::CreateCellCycleModel()
{
    assert(mpCell!=NULL);
    // calls a cheeky version of the constructor which makes the new cell 
    // cycle model the same as the old one - not a dividing copy at this time.
    // unless the parent cell has just divided.        
    return new OxygenBasedCellCycleModel(mpOdeSystem, mpCell->GetMutationState(), 
                                         mBirthTime, mLastTime, mReadyToDivide, mDivideTime);
}


void OxygenBasedCellCycleModel::Initialise()
{
    assert(mpOdeSystem==NULL);
    assert(mpCell!=NULL);    
    mpOdeSystem = new Alarcon2004OxygenBasedCellCycleOdeSystem(CellwiseData<2>::Instance()->GetValue(mpCell,0), mpCell->GetMutationState());
    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());  
}    
 
