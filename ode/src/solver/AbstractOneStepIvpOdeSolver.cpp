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
    
    int numberOfTimeSamples;
    double currentTime;
    
    numberOfTimeSamples = 0;
    
    currentTime = startTime;
    
    while (currentTime < endTime)
    {
        numberOfTimeSamples++;
        
        if (startTime+numberOfTimeSamples*timeSampling >= endTime)
        {
            currentTime = endTime;
        }
        else
        {
            currentTime = startTime+numberOfTimeSamples*timeSampling;
        }
    }
    
    // setup solutions if output is required
    
    OdeSolution solutions;
    
    solutions.SetNumberOfTimeSteps(numberOfTimeSamples);
    solutions.rGetSolutions().push_back(rYValues);
    solutions.rGetTimes().push_back(startTime);
    
    // Solve the ODE system
    
    int timeStepNumber = 0;
    
    currentTime = startTime;
    
    double toTime;
    
    while (currentTime < endTime)
    {
        timeStepNumber++;
        
        toTime = startTime+timeStepNumber*timeSampling;
        
        if (toTime >= endTime)
        {
            toTime = endTime;
        }
        
        Solve(pAbstractOdeSystem, rYValues, currentTime, toTime, timeStep);
        
        currentTime = toTime;
        if ( CalculateStoppingEvent(pAbstractOdeSystem, rYValues, currentTime) == true )
        {
            endTime=currentTime;
            solutions.SetNumberOfTimeSteps(timeStepNumber);
        }
        
        // write current solution into solutions
        solutions.rGetSolutions().push_back(rYValues);
        // Push back new time into the time solution vector
        solutions.rGetTimes().push_back(currentTime);
        
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
    
    int timeStepNumber = 0;
    
    double currentTime = startTime;
    
    while (currentTime < endTime)
    {
        timeStepNumber++;
        
        // Determine what the value time step should really be like
        
        if (startTime+timeStepNumber*timeStep >= endTime)
        {
            realTimeStep = endTime-currentTime;
        }
        
        // Function that calls the appropriate one-step solver
//        std::cout << rYValues[0] << "\n";
//        std::cout << pAbstractOdeSystem->GetNumberOfStateVariables()<<"\n";
//        std::cout.flush();
        rYValues = CalculateNextYValue(pAbstractOdeSystem,
                                       realTimeStep,
                                       currentTime,
                                       rYValues);
                                       
        // Determine the new current time
        
        if (realTimeStep < timeStep)
        {
            currentTime = endTime;
        }
        else
        {
            currentTime = startTime+timeStepNumber*timeStep;
        }
    }
}
