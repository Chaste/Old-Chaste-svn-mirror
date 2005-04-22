#include "AdamsBashforthIvpOdeSolver.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"

#include <iostream>
#include <vector>
#include <cassert>

/*
 * Adams-Bashforth Method Initial Value Problem Ordinary Differential Equation Solver
 * 
 * Solves a system of ODEs for a given time range and time step.
 * To be used in the form:
 * 
 * AdamsBashforthIvpOdeSolver mySolver
 * OdeSolution solution=mySolver->Solve(pMyOdeSystem, StartTime, EndTime, TimeStep, yInit);
 * 
 * where:
 * pMyOdeSystem is a pointer to a specific instance of a subclass of AbstractOdeSystem
 * (this defines the derivatives of the system)
 * the times are all doubles
 * yInit is a std::vector of doubles with initial values for all unknowns
 * 
 * OdeSolution is an object containing an integer of the number of equations, 
 * a std::vector of times and a std::vector of std::vectors of the solution 
 * of the ODE system at those times
*/

OdeSolution AdamsBashforthIvpOdeSolver::Solve(AbstractOdeSystem* pAbstractOdeSystem, 
				double startTime,
				double endTime,
				double timeStep,
				std::vector<double> initialConditions)
{

    int num_equations = pAbstractOdeSystem->mNumberOfEquations;
    
    // Assert that the size of Initial Conditions vector = number of equations.
    assert(initialConditions.size()==num_equations);	
    
    // Assert that the timestep does not exceed the time interval.
    assert(timeStep <= endTime - startTime);
    
    // Assert that we actually have a time interval > 0 .
    assert(endTime > startTime);
    
    // Assert that we  have a timestep > 0 .
    assert(timeStep > 0.0);
    
    int num_timesteps = ((int) ((endTime - startTime)/timeStep));   
    double last_timestep = endTime - ((double) num_timesteps)*timeStep;
    
	OdeSolution solutions;
	solutions.mNumberOfTimeSteps = num_timesteps;
		
	solutions.mSolutions.push_back(initialConditions);
	solutions.mTime.push_back(startTime);
	
	std::vector<double> row(num_equations);	
	std::vector<double> dy(num_equations);
	
	row=initialConditions;
	
	std::vector<std::vector<double> > temp;
	
	std::vector<double> k1(num_equations);
	std::vector<double> k2(num_equations);
	std::vector<double> k3(num_equations);
	std::vector<double> k4(num_equations);
	
	std::vector<double> yk2(num_equations);
	std::vector<double> yk3(num_equations);
	std::vector<double> yk4(num_equations);
	
	for(int timeindex=0; timeindex<3; timeindex++)
	{
		// Apply RungeKutta4's method first three timesteps, in order to 
		// maintain fourth order accuracy of Adams-Bashforth method
		
        k1 = pAbstractOdeSystem->EvaluateYDerivatives(solutions.mTime[timeindex],row);
        
        for(int i=0;i<num_equations; i++) 
		{
			yk2[i] = row[i] + 0.5*k1[i];		
		}
        k2 = pAbstractOdeSystem->EvaluateYDerivatives(solutions.mTime[timeindex]+0.5*timeStep,yk2);
		
		for(int i=0;i<num_equations; i++) 
		{
			yk3[i] = row[i] + 0.5*k2[i];		
		}
        k3 = pAbstractOdeSystem->EvaluateYDerivatives(solutions.mTime[timeindex]+0.5*timeStep,yk3);        

		for(int i=0;i<num_equations; i++) 
		{
			yk4[i] = row[i] + k3[i];		
		}
        k4 = pAbstractOdeSystem->EvaluateYDerivatives(solutions.mTime[timeindex]+timeStep,yk4);                
		
		for(int i=0;i<num_equations; i++) 
		{
			row[i] = row[i] + timeStep*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;		
		}
		
		solutions.mSolutions.push_back(row);	
		solutions.mTime.push_back(solutions.mTime[timeindex]+timeStep);
		temp.push_back(dy);
	}
	    
	// Apply Adams-Bashforth method
    for (int timeindex=3; timeindex<num_timesteps; timeindex++)
    {
    	
        
    	dy = pAbstractOdeSystem->EvaluateYDerivatives(solutions.mTime[timeindex],row);
        
  		for(int i=0;i<num_equations; i++) 
		{		
			row[i] = row[i] + (timeStep/24.0)*(55.0*dy[i] - 59.0*temp[timeindex-1][i] + 37.0*temp[timeindex-2][i] - 9.0*temp[timeindex-3][i]);		
		}
		
		solutions.mSolutions.push_back(row);
		solutions.mTime.push_back(solutions.mTime[timeindex]+timeStep);
		temp.push_back(dy);
	}
	
	// Extra step to get to exactly endTime
	int timeindex = num_timesteps;
	if(last_timestep > 0.00001)
	{	
		solutions.mNumberOfTimeSteps=num_timesteps+1;
		dy = pAbstractOdeSystem->EvaluateYDerivatives(solutions.mTime[num_timesteps],row);
		for(int i=0;i<num_equations; i++) 
		{				
			row[i] = row[i] + (timeStep/24.0)*(55.0*dy[i] - 59.0*temp[timeindex-1][i] + 37.0*temp[timeindex-2][i] - 9.0*temp[timeindex-3][i]);		
		}
		
		solutions.mSolutions.push_back(row);
		solutions.mTime.push_back(solutions.mTime[num_timesteps]+last_timestep);
	}
	
			
	return solutions;
}
