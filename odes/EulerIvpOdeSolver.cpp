/**
 * Concrete EulerIvpOdeSolver class. 
*/
#include "EulerIvpOdeSolver.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"

#include <iostream>
#include <vector>
#include <cassert>

/**
 * Solves a system of ODEs using the Forward Euler method
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
 * system at those times
 * 
 * To be used in the form:
 * 
 * EulerIvpOdeSolver mySolver;
 * OdeSolution solution=mySolver->Solve(pMyOdeSystem, StartTime, EndTime, TimeStep, yInit);
 *  
*/

std::vector<double> EulerIvpOdeSolver::CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
										double timeStep,
										double time,
										std::vector<double> currentYValue)
{
	int num_equations = pAbstractOdeSystem->mNumberOfEquations;
	std::vector<double> dy(num_equations);
	std::vector<double> next_y_value(num_equations);
    dy = pAbstractOdeSystem->EvaluateYDerivatives(time, currentYValue);
    
	for(int i=0;i<num_equations; i++) 
	{
		next_y_value[i] = currentYValue[i] + timeStep*dy[i];		
	}
	return next_y_value;
}
