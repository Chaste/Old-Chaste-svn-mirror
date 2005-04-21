#include "EulerIvpOdeSolver.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"

#include <iostream>
#include <vector>

/*
 * OdeSolution is an object containing an integer of the number of equations, 
 * a std::vector of times and a std::vector of std::vectors of the solution 
 * of the ODE system at those times
*/

OdeSolution EulerIvpOdeSolver::Solve(AbstractOdeSystem* pAbstractOdeSystem, 
				double startTime,
				double endTime,
				double timeStep,
				std::vector<double> initialConditions)
{

    int num_equations = pAbstractOdeSystem->mNumberOfEquations;
    int num_timesteps = ((int) ((endTime - startTime)/timeStep));
    
    double last_timestep = endTime - ((double) num_timesteps)*timeStep;
    
	OdeSolution solutions;
	solutions.mNumberOfTimeSteps = num_timesteps;
	// (num_timesteps)(num_equations)
		
	solutions.mSolutions.push_back(initialConditions);
	solutions.mTime.push_back(startTime);
	
	std::vector<double> row(num_equations);	
	std::vector<double> dy(num_equations);
	
	row=initialConditions;
	
	for(int timeindex=0;timeindex<num_timesteps;timeindex++)
	{
		// Apply Euler's method
        dy = pAbstractOdeSystem->EvaluateYDerivatives(solutions.mTime[timeindex],row);
        
		for(int i=0;i<num_equations; i++) 
		{
			row[i] = row[i] + timeStep*dy[i];		
		}
		
		solutions.mSolutions.push_back(row);
		
		solutions.mTime.push_back(solutions.mTime[timeindex]+timeStep);
	}
	
	// Extra step to get to exactly endTime
	if(last_timestep>0.00001)
	{	
		solutions.mNumberOfTimeSteps=num_timesteps+1;
		dy = pAbstractOdeSystem->EvaluateYDerivatives(solutions.mTime[num_timesteps+1],row);
		for(int i=0;i<num_equations; i++) 
		{
			row[i] = row[i] + last_timestep*dy[i];		
		}
		solutions.mSolutions.push_back(row);
		
		solutions.mTime.push_back(solutions.mTime[num_timesteps]+last_timestep);
	
	}
	
			
	return solutions;
}
