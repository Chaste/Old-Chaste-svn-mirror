#include "AbstractOneStepIvpOdeSolver.hpp"
#include "TimeStepper.hpp"
#include <cassert>
#include <math.h>


const double smidge=1e-10;

OdeSolution AbstractOneStepIvpOdeSolver::Solve(AbstractOdeSystem* pAbstractOdeSystem,
                                               std::vector<double>& rYValues,
                                               double startTime,
                                               double endTime,
                                               double timeStep,
                                               double timeSampling)
{
    assert(rYValues.size()==pAbstractOdeSystem->GetNumberOfStateVariables());
    assert(endTime > startTime);
    assert(timeStep > 0.0);
    assert(timeSampling >= timeStep);
    
    mStoppingEventOccured = false;
    if ( pAbstractOdeSystem->CalculateStoppingEvent(startTime, rYValues) == true )
    {
        EXCEPTION("Stopping event is true for initial condition");
    } 
    // Determine the number of time steps that will be required to solve the
    // ODE system (note that the current algorithm accounts for any potential
    // floating point error)
    unsigned guess_number_of_time_samples = (unsigned) ceil((endTime - startTime)/timeSampling);

    // setup solutions if output is required
    OdeSolution solutions;
    //Set number of time steps will be duplicated below
    solutions.SetNumberOfTimeSteps(guess_number_of_time_samples);
    solutions.rGetSolutions().push_back(rYValues);
    solutions.rGetTimes().push_back(startTime);
    
    // Allocate working memory
    std::vector<double> working_memory(rYValues.size());
    
    // Solve the ODE system
    unsigned time_step_number = 0;
    TimeStepper stepper(startTime, endTime, timeSampling);
    while ( !stepper.IsTimeAtEnd() && !mStoppingEventOccured )
    {
        InternalSolve(pAbstractOdeSystem, rYValues, working_memory, stepper.GetTime(), stepper.GetNextTime(), timeStep);
        stepper.AdvanceOneTimeStep();
        time_step_number++;
        // write current solution into solutions
        solutions.rGetSolutions().push_back(rYValues);
        // Push back new time into the time solution vector
        if ( mStoppingEventOccured )
        {
            solutions.rGetTimes().push_back(mStoppingTime);
        }
        else
        {
            solutions.rGetTimes().push_back(stepper.GetTime());
        }
    }
    
    //This line is here because the above first "reservation" time loop
    //may have different behaviour under optimisation than the second
    //"calculation" time loop.
    solutions.SetNumberOfTimeSteps(time_step_number);
    return solutions;
}

void AbstractOneStepIvpOdeSolver::Solve(AbstractOdeSystem* pAbstractOdeSystem,
                                        std::vector<double>& rYValues,
                                        double startTime,
                                        double endTime,
                                        double timeStep)
{
    // Allocate working memory
    std::vector<double> working_memory(rYValues.size());
    // And solve...
    InternalSolve(pAbstractOdeSystem, rYValues, working_memory, startTime, endTime, timeStep);
}


void AbstractOneStepIvpOdeSolver::InternalSolve(AbstractOdeSystem* pAbstractOdeSystem,
                                                std::vector<double>& rYValues,
                                                std::vector<double>& rWorkingMemory,
                                                double startTime,
                                                double endTime,
                                                double timeStep)
{
    TimeStepper stepper(startTime, endTime, timeStep);
    // Solve the ODE system
    
    // Which of our vectors holds the current solution?
    // If this is true, it's in rYValues, otherwise it's in rWorkingMemory.
    bool curr_is_curr = false;
    
    // should never get here if this bool has been set to true;
    assert(!mStoppingEventOccured);
    while ( !stepper.IsTimeAtEnd() && !mStoppingEventOccured )
    {
        curr_is_curr = not curr_is_curr;
        // Function that calls the appropriate one-step solver
        CalculateNextYValue(pAbstractOdeSystem,
                            stepper.GetNextTimeStep(),
                            stepper.GetTime(),
                            curr_is_curr ? rYValues : rWorkingMemory,
                            curr_is_curr ? rWorkingMemory : rYValues);
        stepper.AdvanceOneTimeStep();
        if ( pAbstractOdeSystem->CalculateStoppingEvent(stepper.GetTime(),
                                                        curr_is_curr ? rWorkingMemory : rYValues) == true )
        {
            mStoppingTime = stepper.GetTime();
            mStoppingEventOccured = true;
        }
    } 
    // Final answer must be in rYValues
    if (curr_is_curr)
    {
        rYValues.assign(rWorkingMemory.begin(), rWorkingMemory.end());
    }
}
