#include "OxygenBasedCellCycleModel.hpp"
#include "Exception.hpp"
#include <iostream>
#include <cassert>
#include <cfloat>

RungeKutta4IvpOdeSolver OxygenBasedCellCycleModel::msSolver;

/**
 * Model constructor.
 *
 * Note that since we provide our own constructor, the compiler will *not*
 * generate a default one for us.
 *
 */

OxygenBasedCellCycleModel::OxygenBasedCellCycleModel()
        : AbstractCellCycleModel(),
          mpOdeSystem(NULL)
{
    SimulationTime* p_sim_time = SimulationTime::Instance();
    if (p_sim_time->IsStartTimeSetUp()==false)
    {
        EXCEPTION("OxygenBasedCellCycleModel is being created but SimulationTime has not been set up");
    }
    mBirthTime = p_sim_time->GetDimensionalisedTime();
    mLastTime = mBirthTime;
    mReadyToDivide = false;
    mDivideTime = DBL_MAX;

}

/**
 * A private constructor for daughter cells called only by the CreateCellCycleModel function
 *
 * @param parentProteinConcentrations a std::vector of doubles of the protein concentrations (see OxygenBasedCellCycleModel)
 * @param birthTime the simulation time when the cell divided (birth time of parent cell)
 */

OxygenBasedCellCycleModel::OxygenBasedCellCycleModel(Alarcon2004OxygenBasedCellCycleOdeSystem* pParentOdeSystem, 
                              const bool& rIsCancerCell, double birthTime, 
                              double lastTime, bool readyToDivide, double divideTime)
        : AbstractCellCycleModel()
{
    if (pParentOdeSystem !=NULL)
    {
        std::vector<double> parent_protein_concs = pParentOdeSystem->rGetStateVariables();
        mpOdeSystem = new Alarcon2004OxygenBasedCellCycleOdeSystem(parent_protein_concs[5], rIsCancerCell);
        
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
    mLastTime = lastTime;
    mReadyToDivide = readyToDivide;
    mDivideTime = divideTime;
}

/**
 * A private constructor for archiving
 *
 * @param parentProteinConcentrations a std::vector of doubles of the protein concentrations (see OxygenBasedCellCycleModel)
 * @param birthTime the simulation time when the cell divided (birth time of parent cell)
 */

OxygenBasedCellCycleModel::OxygenBasedCellCycleModel(const std::vector<double>& rParentProteinConcentrations, 
                              const bool& rIsCancerCell, double birthTime, 
                              double lastTime, bool readyToDivide, double divideTime)
        : AbstractCellCycleModel()
{
    mpOdeSystem = new Alarcon2004OxygenBasedCellCycleOdeSystem(rParentProteinConcentrations[5], rIsCancerCell);
    
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
    mLastTime = lastTime;
    mReadyToDivide = readyToDivide;
    mDivideTime = divideTime;
}


OxygenBasedCellCycleModel::~OxygenBasedCellCycleModel()
{
    delete mpOdeSystem;
}

/**
 * Resets the oxygen-based model to the start of the cell cycle (this model does not cycle naturally)
 * Cells are given a new birth time and cell cycle proteins are reset.
 * Note that the wnt pathway proteins maintain their current values.
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
            
    // Once MeinekeCryptCell has been re-factored, the new daughter TumourCell class
    // can have a GetIsCancerCell() method. Hard-code it for the time being.    
    mpOdeSystem->SetIsCancerCell(false);
    
    double current_time = SimulationTime::Instance()->GetDimensionalisedTime();
    
    if (current_time>mLastTime)
    {
        // feed this time step's oxygen concentration into the solver as a constant over this timestep.
        double dt = 0.0001; // Needs to be this precise to stop crazy errors whilst we are still using rk4.
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
 * Returns the protein concentrations at the current time (useful for tests)
 *
 * NB: Will copy the vector - you can't use this to modify the concentrations.
 */
std::vector<double>  OxygenBasedCellCycleModel::GetProteinConcentrations() const
{
    assert(mpOdeSystem!=NULL);
    return mpOdeSystem->rGetStateVariables();
}

/**
 * Returns a new WntCellCycleModel created with the correct initial conditions.
 *
 * Should be called just after the parent cell cycle model has been .Reset().
 */
AbstractCellCycleModel* OxygenBasedCellCycleModel::CreateCellCycleModel()
{
    assert(mpCell!=NULL);
    // calls a cheeky version of the constructor which makes the new cell 
    // cycle model the same as the old one - not a dividing copy at this time.
    // unless the parent cell has just divided.
    
    // Once MeinekeCryptCell has been re-factored, the new daughter TumourCell class
    // can have a GetIsCancerCell() method. Hard-code it for the time being.
    return new OxygenBasedCellCycleModel(mpOdeSystem, 
                                         false, mBirthTime, mLastTime, 
                                         mReadyToDivide, mDivideTime);
}

/**
 * Returns a new OxygenBasedCellCycleModel created with the correct initial conditions.
 *
 * Should only be used in tests
 *
 * @param birthTime the simulation time when the cell was born
 */
void OxygenBasedCellCycleModel::SetBirthTime(double birthTime)
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
void OxygenBasedCellCycleModel::SetProteinConcentrationsForTestsOnly(double lastTime, std::vector<double> proteinConcentrations)
{
    assert(mpOdeSystem!=NULL);
    assert(proteinConcentrations.size()==mpOdeSystem->rGetStateVariables().size());
    mLastTime = lastTime;
    mpOdeSystem->SetStateVariables(proteinConcentrations);
}

void OxygenBasedCellCycleModel::Initialise()
{
    assert(mpOdeSystem==NULL);
    assert(mpCell!=NULL);
    
    double oxygen_concentration = CellwiseData<2>::Instance()->GetValue(mpCell,0);
    // Once MeinekeCryptCell has been re-factored, the new daughter TumourCell class
    // can have a GetIsCancerCell() method. Hard-code it for the time being.
    mpOdeSystem = new Alarcon2004OxygenBasedCellCycleOdeSystem(oxygen_concentration, false);
    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());  
}    


