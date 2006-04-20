/**
 * Concrete AdamsBashforthIvpOdeSolver class. Sub-class of AbstractIvpOdeSolver.hpp
*/
#include "AdamsBashforthIvpOdeSolver.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"
#include "Exception.hpp"

//#include <iostream>
#include <vector>
#include <cassert>

/**
 * Solves a system of ODEs using the Adams-Bashforth Method Initial Value Problem Ordinary Differential Equation Solver
 * 
 * @param pAbstractOdeSystem points to the concrete ODE system to be solved
 * @param rYValues a standard vector specifying the intial condition 
 * of each solution variable in the system 
 * @param startTime the time at which the initial conditions are specified
 * @param endTime the time to which the system should be solved and the solution 
 * returned
 * @param timeStep the time interval to be used by the solver
 * @param timeSampling the time interval for generating the solution
 * 
 * 
 * @return OdeSolution is an object containing an integer of the number of 
 * equations, a std::vector of times and a std::vector of std::vectors where 
 * each of those vectors contains the solution for one variable of the ODE 
 * system at those times
 * 
 * To be used in the form:
 * 
 * AdamsBashforthIvpOdeSolver solver
 * OdeSolution solution=solver.Solve(pAbstractOdeSystem, rYValues, startTime, endTime, timeStep, timeStep);
 *  
 * Note: at this point in time, timeSampling should be EXACTLY the same
 *       as timeStep or it will NOT work! This should not be too much of
 *       an issue since we don't really make use of this integrator.
 *       However, should we need it, we would then need to modify this
 *       class to account for timeSampling.
 *  
 */

OdeSolution AdamsBashforthIvpOdeSolver::Solve(AbstractOdeSystem* pAbstractOdeSystem, 
                                              std::vector<double>& rYValues,
				                              double startTime,
                              				  double endTime,
				                              double timeStep,
                                              double timeSampling)
{
    unsigned num_equations = pAbstractOdeSystem->GetNumberOfStateVariables();
    
    // Assert that the size of Initial Conditions vector = number of equations.
    assert(rYValues.size()==num_equations);	
    
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

    // Determine the number of time steps and make sure that we have at least 4 of them

    int numberOfTimeSteps = 0;
    
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

    if (numberOfTimeSteps <= 4)
    {
        throw Exception("A multi-step solver needs at least 4 time steps.");
    }

	std::vector<double> dy(num_equations);
	
	std::vector<std::vector<double> > temp;
	
	std::vector<double> dyRK4(num_equations);
	std::vector<double> k1(num_equations);
	std::vector<double> k2(num_equations);
	std::vector<double> k3(num_equations);
	std::vector<double> k4(num_equations);
	
	std::vector<double> yk2(num_equations);
	std::vector<double> yk3(num_equations);
	std::vector<double> yk4(num_equations);
	
	for(unsigned int timeindex=0; timeindex<3; timeindex++)
	{
		// Apply RungeKutta4's method first three timesteps, in order to 
		// maintain fourth order accuracy of Adams-Bashforth method
		
        dy = dyRK4 = pAbstractOdeSystem->EvaluateYDerivatives(solutions.rGetTimes()[timeindex],rYValues);
        
        for(unsigned int i=0;i<num_equations; i++) 
		{
			k1[i] = timeStep*dyRK4[i];
			yk2[i] = rYValues[i] + 0.5*k1[i];		
		}
        dyRK4 = pAbstractOdeSystem->EvaluateYDerivatives(solutions.rGetTimes()[timeindex]+0.5*timeStep,yk2);
		
		for(unsigned int i=0;i<num_equations; i++) 
		{
			k2[i] = timeStep*dyRK4[i];
			yk3[i] = rYValues[i] + 0.5*k2[i];		
		}
        dyRK4 = pAbstractOdeSystem->EvaluateYDerivatives(solutions.rGetTimes()[timeindex]+0.5*timeStep,yk3);        

		for(unsigned int i=0;i<num_equations; i++) 
		{
			k3[i] = timeStep*dyRK4[i];
			yk4[i] = rYValues[i] + k3[i];		
		}
        dyRK4 = pAbstractOdeSystem->EvaluateYDerivatives(solutions.rGetTimes()[timeindex]+timeStep,yk4);                
		
		for(unsigned int i=0;i<num_equations; i++) 
		{
			k4[i] = timeStep*dyRK4[i];
			rYValues[i] = rYValues[i] + (k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;
			
		}

//!!! Think about whether OdeSolutions need updating, based on the value of timeSampling
		
		solutions.rGetSolutions().push_back(rYValues);	
		solutions.rGetTimes().push_back(solutions.rGetTimes()[timeindex]+timeStep);
		temp.push_back(dy);
	}

    double realTimeStep = timeStep;

    int timeStepNumber = 2;
    
    currentTime = startTime+timeStepNumber*timeStep;

	// Apply Adams-Bashforth method

    while (currentTime < endTime)
    {
        timeStepNumber++;

        // Determine what the value time step should really be like
        
        if (startTime+timeStepNumber*timeStep >= endTime)
        {
            realTimeStep = endTime-currentTime;
        }

        // Function that calls the appropriate one-step solver
        
        dy = pAbstractOdeSystem->EvaluateYDerivatives(solutions.rGetTimes()[timeStepNumber],rYValues);

        for(unsigned int i=0;i<num_equations; i++) 
        {       
            rYValues[i] = rYValues[i] + (realTimeStep/24.0)*(55.0*dy[i] - 59.0*temp[timeStepNumber-1][i] + 37.0*temp[timeStepNumber-2][i] - 9.0*temp[timeStepNumber-3][i]);        
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

        // Update OdeSolution
//!!! This does NOT take into account timeSampling
        solutions.rGetSolutions().push_back(rYValues);
        solutions.rGetTimes().push_back(currentTime);
        temp.push_back(dy);
    }
	
	return solutions;
}

/**
 * This method is required, since it occurs in the abstract base class.
 * However, it only makes much sense for the one step solvers.
 * Hence here it just behaves as the other solve method, and discards the OdeSolution.
 */
void AdamsBashforthIvpOdeSolver::Solve(AbstractOdeSystem* pAbstractOdeSystem, 
                                       std::vector<double>& rYValues,
                                       double startTime,
                                       double endTime,
                                       double timeStep)
{
    OdeSolution solution = Solve(pAbstractOdeSystem, rYValues, startTime, endTime, timeStep, timeStep);
    rYValues = solution.rGetSolutions()[solution.GetNumberOfTimeSteps()];
}
