#include "AbstractOneStepIvpOdeSolver.hpp"
#include <cassert>

OdeSolution AbstractOneStepIvpOdeSolver::Solve(AbstractOdeSystem* pAbstractOdeSystem, 
				              double startTime,
				              double endTime,
				              double timeStep,
				              std::vector<double> initialConditions)

{	
    int num_equations = pAbstractOdeSystem->mNumberOfEquations;
    
    // Assert that the size of initialConditions vector = number of equations.
    assert(initialConditions.size()==num_equations);	
    
    // Assert that the timestep does not exceed the time interval.
    assert(timeStep < endTime - startTime  + 1e-10);
    
    // Assert that we actually have a time interval > 0 .
    assert(endTime > startTime);
    
    // Assert that we  have a timestep > 0 .
    assert(timeStep > 0.0);
    
    int num_timesteps = ((int) ((endTime - startTime)/timeStep));
    
    double last_timestep = endTime - ((double) num_timesteps)*timeStep - startTime;
    assert(last_timestep < timeStep + 1e-10); 
    
	OdeSolution solutions;
	solutions.mNumberOfTimeSteps = num_timesteps;
		
	solutions.mSolutions.push_back(initialConditions);
	solutions.mTime.push_back(startTime);
	
	std::vector<double> row(num_equations);	
	std::vector<double> dy(num_equations);
	
	row=initialConditions;
	
	for(int timeindex=0;timeindex<num_timesteps;timeindex++)
	{
		
		row = CalculateNextYValue(pAbstractOdeSystem,
									timeStep,
									solutions.mTime[timeindex],
									row);
		
		solutions.mSolutions.push_back(row);
		// Push back new time into the time solution vector
		solutions.mTime.push_back(solutions.mTime[timeindex]+timeStep);
	}
	
	// Extra step to get to exactly endTime
	if(last_timestep > (0.000001 * timeStep))
	{	
		solutions.mNumberOfTimeSteps=num_timesteps+1;
  		
  		row = CalculateNextYValue(pAbstractOdeSystem,
  									last_timestep, 
  									solutions.mTime[num_timesteps],
  									row);
  		
		solutions.mSolutions.push_back(row);
		
		solutions.mTime.push_back(solutions.mTime[num_timesteps]+last_timestep);
	
	}
				
	return solutions;
}
