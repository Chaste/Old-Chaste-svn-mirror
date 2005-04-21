#include "EulerIvpOdeSolver.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"

//#include "petscvec.h"
#include <iostream>
#include <vector>

// OdeSolution is an object containing PETSc vectors of time and solution of the ODE system at those times

OdeSolution EulerIvpOdeSolver::Solve(AbstractOdeSystem* pAbstractOdeSystem, 
				double startTime,
				double endTime,
				double timeStep,
				std::vector<double> initialConditions)
{

    int num_equations = pAbstractOdeSystem->mNumberOfEquations;
    int num_timesteps = ((int) ((endTime - startTime)/timeStep));
    
    double last_timestep = endTime - ((double) num_timesteps)*timeStep;
    
	// Convert Initial conditions from PETSc to normal vector
	OdeSolution solutions;
	// (num_timesteps)(num_equations)
	
	std::vector<double> row(num_equations);	
	std::vector<double> dy(num_equations);
	
	row=initialConditions;
	
	solutions.mSolutions.push_back(initialConditions);

	
	solutions.mTime.push_back(startTime);
	
	for(int timeindex=0;timeindex<num_timesteps;timeindex++)
	{
		// Work out next time step
		// EULER METHOD HERE ************************************		
        dy = pAbstractOdeSystem->EvaluateYDerivatives(solutions.mTime[timeindex],row);
		for(int i=0;i<num_equations; i++) 
		{
			row[i] = row[i] + timeStep*dy[i];		
		}
		solutions.mSolutions.push_back(row);
		
		solutions.mTime.push_back(solutions.mTime[timeindex]+timeStep);
	}
	
	
	

	// Create the PETSc vector of times		
		solutions.mNumberOfTimeSteps = num_timesteps;
		return solutions;
}
