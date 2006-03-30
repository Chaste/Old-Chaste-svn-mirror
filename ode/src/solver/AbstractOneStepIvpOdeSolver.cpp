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
    // setup solutions if output is required
    if(!bUseStateVariable)
    {
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

	
	for(int time_index=0;time_index<num_timesteps;time_index++)
	{
		// Function that calls the appropriate one-step solver
		//row = CalculateNextYValue(pAbstractOdeSystem,
				//					timeStep,
				//					startTime + time_index*timeStep, //solutions.mTime[timeindex],
				//					row);
		
        if(!bUseStateVariable)
        {
            row = CalculateNextYValue(pAbstractOdeSystem,
                                    timeStep,
                                    startTime + time_index*timeStep, //solutions.mTime[timeindex],
                                    row);
            // write current solution into solutions
		    solutions.mSolutions.push_back(row);
		    // Push back new time into the time solution vector
            solutions.mTime.push_back(startTime+(time_index+1)*timeStep);
        }
        else
        {
            CalculateNextYValue(pAbstractOdeSystem,
                                    timeStep,
                                    startTime + time_index*timeStep);
            //pAbstractOdeSystem->SetStateVariables(row);
        }
	}


	// Extra step to get to exactly endTime
	if(last_timestep > (0.000001 * timeStep))
	{	
  		//row = CalculateNextYValue(pAbstractOdeSystem,
  		//							last_timestep, 
  		//							startTime+num_timesteps*timeStep, //solutions.mTime[num_timesteps],
  		//							row);


  		if(!bUseStateVariable)
        {
            row = CalculateNextYValue(pAbstractOdeSystem,
                                    last_timestep, 
                                    startTime+num_timesteps*timeStep, //solutions.mTime[num_timesteps],
                                    row);
            
		    solutions.mSolutions.push_back(row);
            solutions.mTime.push_back(solutions.mTime[num_timesteps]+last_timestep);
        }
        else
        {
            
            CalculateNextYValue(pAbstractOdeSystem,
                                    last_timestep, 
                                    startTime+num_timesteps*timeStep);
            //pAbstractOdeSystem->SetStateVariables(row);
        }
	}
	
    			
	return solutions;
}
