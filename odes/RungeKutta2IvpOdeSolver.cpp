/**
 * Concrete RungeKutta2IvpOdeSolver class. 
*/
#include "RungeKutta2IvpOdeSolver.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"

#include <iostream>
#include <vector>
#include <cassert>

/**
 * Solves a system of ODEs using the Runge Kutta 2nd Order Initial Value Problem Ordinary Differential Equation Solver
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
 * RungeKutta2IvpOdeSolver mySolver
 * OdeSolution solution=mySolver->Solve(pMyOdeSystem, StartTime, EndTime, TimeStep, yInit);
 *  
*/

std::vector<double> RungeKutta2IvpOdeSolver::CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                                        double timeStep,
                                        double time,
                                        std::vector<double> currentYValue)
{

    int num_equations = pAbstractOdeSystem->mNumberOfEquations;
    
    
		// Apply Runge-Kutta 2nd Order method
        
        // Work out k1
        std::vector<double> k1(num_equations);
        std::vector<double> dy(num_equations);
        std::vector<double> next_y_value(num_equations);
        dy = pAbstractOdeSystem->EvaluateYDerivatives(time,currentYValue);
		for(int i=0;i<num_equations; i++) 
		{
			k1[i] = timeStep*dy[i];
			k1[i] = k1[i]/2.0+currentYValue[i];
		}
		// Work out k2				
		std::vector<double> k2(num_equations);
		dy = pAbstractOdeSystem->EvaluateYDerivatives(time+timeStep/2.0,k1);
		for(int i=0;i<num_equations; i++) 
		{
			k2[i] = timeStep*dy[i];
			next_y_value[i]=currentYValue[i]+k2[i];
		}
		
		return next_y_value;
}
