#include "AbstractOneStepIvpOdeSolver.hpp"
#include <cassert>
#include <iostream>


OdeSolution AbstractOneStepIvpOdeSolver::Solve(AbstractOdeSystem* pAbstractOdeSystem,
                                               std::vector<double>& rYValues,
                                               double startTime,
                                               double endTime,
                                               double timeStep,
                                               double timeSampling)
{
    // assert the size of the rYValues vector is correct
    assert(rYValues.size()==pAbstractOdeSystem->GetNumberOfStateVariables());
    
    // Assert that we actually have a time interval > 0 .
    assert(endTime > startTime);
    
    // Assert that we  have a timestep > 0 .
    assert(timeStep > 0.0);
    
    // Assert that we  have a timesampling > 0 and >= timestep
    assert(timeSampling >= timeStep);
    
    // Determine the number of time steps that will be required to solve the
    // ODE system (note that the current algorithm accounts for any potential
    // floating point error)
    
    mStoppingEventOccured = false;
    if ( pAbstractOdeSystem->CalculateStoppingEvent(startTime, rYValues) == true )
    {
        EXCEPTION("Stopping event is true for initial condition");
    }
 
    int number_of_time_samples;
    double current_time;
    
    number_of_time_samples = 0;
    
    current_time = startTime;
    
    while (current_time < endTime)
    {
        number_of_time_samples++;
        
        if (startTime+number_of_time_samples*timeSampling >= endTime)
        {
            current_time = endTime;
        }
        else
        {
            current_time = startTime+number_of_time_samples*timeSampling;
        }
    }
    
    // setup solutions if output is required
    
    OdeSolution solutions;
    
    solutions.SetNumberOfTimeSteps(number_of_time_samples);
    solutions.rGetSolutions().push_back(rYValues);
    solutions.rGetTimes().push_back(startTime);
    
    // Solve the ODE system
    
    int time_step_number = 0;
    
    current_time = startTime;
    
    double to_time;
    

    while ( (current_time < endTime) && (!mStoppingEventOccured) )
    {
        time_step_number++;
        
        to_time = startTime+time_step_number*timeSampling;
        
        if (to_time >= endTime)
        {
            to_time = endTime;
        }
        
        Solve(pAbstractOdeSystem, rYValues, current_time, to_time, timeStep);
        
        current_time = to_time;
        if ( mStoppingEventOccured == true )
        {
            current_time = mStoppingTime;
            endTime = current_time;
            solutions.SetNumberOfTimeSteps(time_step_number);
        }
        
        // write current solution into solutions
        solutions.rGetSolutions().push_back(rYValues);
        // Push back new time into the time solution vector
        solutions.rGetTimes().push_back(current_time);
    }
    
    return solutions;
}

void AbstractOneStepIvpOdeSolver::Solve(AbstractOdeSystem* pAbstractOdeSystem,
                                        std::vector<double>& rYValues,
                                        double startTime,
                                        double endTime,
                                        double timeStep)
{
    // Solve the ODE system
    
    double realTimeStep = timeStep;
    
    int time_step_number = 0;
    
    double current_time = startTime;
        
    //should never get here if this bool has been set to true;
    assert(!mStoppingEventOccured);
    
    while( (current_time < endTime) && (!mStoppingEventOccured) )
    {
        time_step_number++;
        
        // Determine what the value time step should really be like
        
        if (startTime+time_step_number*timeStep >= endTime)
        {
            realTimeStep = endTime-current_time;
        }
        
        // Function that calls the appropriate one-step solver
//        std::cout << rYValues[0] << "\n";
//        std::cout << pAbstractOdeSystem->GetNumberOfStateVariables()<<"\n";
//        std::cout.flush();
        rYValues = CalculateNextYValue(pAbstractOdeSystem,
                                       realTimeStep,
                                       current_time,
                                       rYValues);
                                       
        // Determine the new current time
        
        if (realTimeStep < timeStep)
        {
            current_time = endTime;
        }
        else
        {
            current_time = startTime+time_step_number*timeStep;
        }
        
        if ( pAbstractOdeSystem->CalculateStoppingEvent(current_time, rYValues) == true )
        {
            mStoppingTime = current_time;
            mStoppingEventOccured = true;
        }
    }
}
