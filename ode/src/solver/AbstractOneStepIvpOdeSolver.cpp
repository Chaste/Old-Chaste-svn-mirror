#include "AbstractOneStepIvpOdeSolver.hpp"
#include "TimeStepper.hpp"
#include <cassert>
#include <math.h>


const double smidge=1e-10;

OdeSolution AbstractOneStepIvpOdeSolver::Solve(AbstractOdeSystem* pOdeSystem,
                                               std::vector<double>& rYValues,
                                               double startTime,
                                               double endTime,
                                               double timeStep,
                                               double timeSampling)
{
    assert(rYValues.size()==pOdeSystem->GetNumberOfStateVariables());
    assert(endTime > startTime);
    assert(timeStep > 0.0);
    assert(timeSampling >= timeStep);
    
    mStoppingEventOccured = false;
    if ( pOdeSystem->CalculateStoppingEvent(startTime, rYValues) == true )
    {
        EXCEPTION("(Solve with sampling) Stopping event is true for initial condition");
    } 
    TimeStepper stepper(startTime, endTime, timeSampling);

    // setup solutions if output is required
    OdeSolution solutions;
    solutions.SetNumberOfTimeSteps(stepper.EstimateTimeSteps());
    solutions.rGetSolutions().push_back(rYValues);
    solutions.rGetTimes().push_back(startTime);
    
    mWorkingMemory.resize(rYValues.size());
    
    // Solve the ODE system
    while ( !stepper.IsTimeAtEnd() && !mStoppingEventOccured )
    {
        InternalSolve(pOdeSystem, rYValues, mWorkingMemory, stepper.GetTime(), stepper.GetNextTime(), timeStep);
        stepper.AdvanceOneTimeStep();
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
    
    // stepper.EstimateTimeSteps may have been an overestimate...
    solutions.SetNumberOfTimeSteps(stepper.GetTimeStepsElapsed());
    return solutions;
}

void AbstractOneStepIvpOdeSolver::Solve(AbstractOdeSystem* pOdeSystem,
                                        std::vector<double>& rYValues,
                                        double startTime,
                                        double endTime,
                                        double timeStep)
{
    assert(rYValues.size()==pOdeSystem->GetNumberOfStateVariables());
    assert(endTime > startTime);
    assert(timeStep > 0.0);
    
    mStoppingEventOccured = false;
    if ( pOdeSystem->CalculateStoppingEvent(startTime, rYValues) == true )
    {
        EXCEPTION("(Solve without sampling) Stopping event is true for initial condition");
    }
    
    // Perhaps resize working memory
    mWorkingMemory.resize(rYValues.size());
    // And solve...
    InternalSolve(pOdeSystem, rYValues, mWorkingMemory, startTime, endTime, timeStep);
}


void AbstractOneStepIvpOdeSolver::InternalSolve(AbstractOdeSystem* pOdeSystem,
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
        CalculateNextYValue(pOdeSystem,
                            stepper.GetNextTimeStep(),
                            stepper.GetTime(),
                            curr_is_curr ? rYValues : rWorkingMemory,
                            curr_is_curr ? rWorkingMemory : rYValues);
        stepper.AdvanceOneTimeStep();
        if ( pOdeSystem->CalculateStoppingEvent(stepper.GetTime(),
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
