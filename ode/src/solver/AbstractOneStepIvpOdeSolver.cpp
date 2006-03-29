#include "AbstractOneStepIvpOdeSolver.hpp"
#include <cassert>

/*
 * Solves a system of ODEs using a specified one-step ODE solver
 * 
 * @param pAbstractOdeSystem points to the concrete ODE system to be solved
 * @param startTime the time at which the initial conditions are specified
 * @param endTime the time to which the system should be solved and the solution 
 * returned
 * @param timeStep the time interval to be used by the solver
 * @param initialConditions a standard vector specifying the intial condition 
 * of each solution variable in the system 
 * 
 * 
 * @return OdeSolution is an object containing an integer of the number of 
 * equations, a std::vector of times and a std::vector of std::vectors where 
 * each of those vectors contains the solution for one variable of the ODE 
 * system at those times.
 */
 
OdeSolution AbstractOneStepIvpOdeSolver::Solve(AbstractOdeSystem* pAbstractOdeSystem, 
				              double startTime,
				              double endTime,
				              double timeStep,
				              std::vector<double> initialConditions)

{	
    // Assert that the size of initialConditions vector = number of equations.
    assert(initialConditions.size()==pAbstractOdeSystem->GetNumberOfStateVariables());	
    
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
	if (last_timestep > (0.000001 * timeStep))
	{
	    // We'll use an extra time step
	    solutions.SetNumberOfTimeSteps(num_timesteps+1);
	}
	else
	{
	    solutions.SetNumberOfTimeSteps(num_timesteps);
	}
		
	solutions.mSolutions.push_back(initialConditions);
	solutions.mTime.push_back(startTime);
	
	std::vector<double> row; // A vector of current Y values.
	
	row=initialConditions;
	
	for(int timeindex=0;timeindex<num_timesteps;timeindex++)
	{
		// Function that calls the appropriate one-step solver
		row = CalculateNextYValue(pAbstractOdeSystem,
									timeStep,
									solutions.mTime[timeindex],
									row);
		
		solutions.mSolutions.push_back(row);
		// Push back new time into the time solution vector
//		solutions.mTime.push_back(solutions.mTime[timeindex]+timeStep);
        solutions.mTime.push_back(startTime+(timeindex+1)*timeStep);
	}
	
	// Extra step to get to exactly endTime
	if(last_timestep > (0.000001 * timeStep))
	{	
	    //solutions.SetNumberOfTimeSteps(num_timesteps+1);
  		
  		row = CalculateNextYValue(pAbstractOdeSystem,
  									last_timestep, 
  									solutions.mTime[num_timesteps],
  									row);
  		
		solutions.mSolutions.push_back(row);
		
		solutions.mTime.push_back(solutions.mTime[num_timesteps]+last_timestep);
	
	}
				
	return solutions;
}

/*
 * Solves a system of ODEs using a specified one-step ODE solver
 * 
 * @param pAbstractOdeSystem points to the concrete ODE system to be solved
 * @param startTime the time at which the initial conditions are specified
 * @param endTime the time to which the system should be solved and the solution 
 * returned
 * @param timeStep the time interval to be used by the solver
 *
 * This version of solve modifies the StateVariables member of pAbstractOdeSystem
 * instead of returning an OdeSolution set. The StateVariables should be initialised 
 * as the initial conditions. 
 */
 
void AbstractOneStepIvpOdeSolver::Solve(AbstractOdeSystem* pAbstractOdeSystem, 
                              double startTime,
                              double endTime,
                              double timeStep)

{   
 
    // Assert that the timestep does not exceed the time interval.
    assert(timeStep < endTime - startTime  + 1e-10);
    
    // Assert that we actually have a time interval > 0 .
    assert(endTime > startTime);
    
    // Assert that we  have a timestep > 0 .
    assert(timeStep > 0.0);
    
    int num_timesteps = ((int) ((endTime - startTime)/timeStep));
    
    double last_timestep = endTime - ((double) num_timesteps)*timeStep - startTime;
    assert(last_timestep < timeStep + 1e-10); 
    
//    OdeSolution solutions;
//    if (last_timestep > (0.000001 * timeStep))
//    {
//        // We'll use an extra time step
//        solutions.SetNumberOfTimeSteps(num_timesteps+1);
//    }
//    else
//    {
//        solutions.SetNumberOfTimeSteps(num_timesteps);
//    }
//        
//    solutions.mSolutions.push_back(initialConditions);
//    solutions.mTime.push_back(startTime);
    
    std::vector<double> row; // A vector of current Y values.
    
    row=pAbstractOdeSystem->GetStateVariables();
    
    for(int time_index=0;time_index<num_timesteps;time_index++)
    {
        // Function that calls the appropriate one-step solver
        row = CalculateNextYValue(pAbstractOdeSystem,
                                    timeStep,
                                    startTime+time_index*timeStep,
                                    row);
        
       pAbstractOdeSystem->SetStateVariables(row);
       
       // solutions.mSolutions.push_back(row);
        // Push back new time into the time solution vector
//      solutsolutions.mTime[num_timesteps],startTime+num_timesteps*timeStepions.mTime.push_back(solutions.mTime[timeindex]+timeStep);
       //solutions.mTime.push_back(startTime+(timeindex+1)*timeStep);
    }
    
    // Extra step to get to exactly endTime
    if(last_timestep > (0.000001 * timeStep))
    {   
        //solutions.SetNumberOfTimeSteps(num_timesteps+1);
        
        row = CalculateNextYValue(pAbstractOdeSystem,
                                    last_timestep, 
                                    endTime,
                                    row);
                                    
        pAbstractOdeSystem->SetStateVariables(row);
        //solutions.mSolutions.push_back(row);
       // solutions.mTime[num_timesteps],startTime+num_timesteps*timeStep
        //solutions.mTime.push_back(solutions.mTime[num_timesteps]+last_timestep);
    
    }
                
    //return solutions;
}
