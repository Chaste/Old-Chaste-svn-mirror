#include "AbstractOneStepIvpOdeSolver.hpp"
#include <cassert>
#include <iostream>


OdeSolution AbstractOneStepIvpOdeSolver::Solve(AbstractOdeSystem* pAbstractOdeSystem, 
                              double startTime,
                              double endTime,
                              double timeStep,
                              std::vector<double> initialConditions /*default is empty vector*/)

{   
    // if no initial conditions, this method just updates the 
    // state variable. Set a boolean here
    bool bUseStateVariable = initialConditions.empty();
    
    if(!bUseStateVariable)
    {    // Assert that the size of initialConditions vector = number of equations.
        assert(initialConditions.size()==pAbstractOdeSystem->GetNumberOfStateVariables());  
    }

    // Assert that we actually have a time interval > 0 .
    assert(endTime > startTime);
    
    // Assert that we  have a timestep > 0 .
    assert(timeStep > 0.0);
    
    int numberOfTimeSteps;
    double currentTime;

    // Determine the number of time steps that will be required to solve the
    // ODE system (note that the current algorithm accounts for any potential
    // floating point error)

    numberOfTimeSteps = 0;
    
    currentTime = startTime;

    while (currentTime < endTime)
    {
        numberOfTimeSteps++;
        
        if (startTime+numberOfTimeSteps*timeStep >= endTime)
        {
            currentTime = endTime;
        }
        else
        {
            currentTime = startTime+numberOfTimeSteps*timeStep;
        }
    }

    // setup solutions if output is required

    OdeSolution solutions;

    if(!bUseStateVariable)
    {
        solutions.SetNumberOfTimeSteps(numberOfTimeSteps);
        solutions.mSolutions.push_back(initialConditions);
        solutions.mTime.push_back(startTime);
    }

    std::vector<double> row; // A vector of current Y values.

    if(!bUseStateVariable)
    {
        // use the input parameter
        row = initialConditions;
    } 
    else
    {
        // use the state variable as the initial condition
        row = pAbstractOdeSystem->GetStateVariables();
    }
    
    // Solve the ODE system

    double realTimeStep = timeStep;

    int timeStepNumber = 0;
    
    currentTime = startTime;

    while (currentTime < endTime)
    {
        timeStepNumber++;

        // Determine what the value time step should really be like
        
        if (startTime+timeStepNumber*timeStep >= endTime)
        {
            realTimeStep = endTime-currentTime;
        }

        // Function that calls the appropriate one-step solver
        
        if(!bUseStateVariable)
        {
            row = CalculateNextYValue(pAbstractOdeSystem,
                                    realTimeStep,
                                    currentTime,
                                    row);
        }
        else
        {
            CalculateNextYValue(pAbstractOdeSystem,
                                    timeStep,
                                    currentTime);
            //pAbstractOdeSystem->SetStateVariables(row);
        }

        // Determine the new current time
    
        if (realTimeStep < timeStep)
        {
            currentTime = endTime;
        }
        else
        {
            currentTime = startTime+timeStepNumber*timeStep;
        }
        
        if (!bUseStateVariable)
        {
            // write current solution into solutions
            solutions.mSolutions.push_back(row);
            // Push back new time into the time solution vector
            solutions.mTime.push_back(currentTime);
        }
    }

    return solutions;
}
