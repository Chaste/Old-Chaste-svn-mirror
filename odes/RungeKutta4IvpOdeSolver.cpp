/**
 * Concrete RungeKutta4IvpOdeSolver class. 
*/
#include "RungeKutta4IvpOdeSolver.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"

#include <iostream>
#include <vector>
#include <cassert>

/**
 * Solves a system of ODEs using the Runge Kutta 4th Order Initial Value Problem Ordinary Differential Equation Solver
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
 * RungeKutta4IvpOdeSolver mySolver
 * OdeSolution solution=mySolver->Solve(pMyOdeSystem, StartTime, EndTime, TimeStep, yInit);
 *  
*/

std::vector<double> RungeKutta4IvpOdeSolver::CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                                        double timeStep,
                                        double time,
                                        std::vector<double> currentYValue)
{
    int num_equations = pAbstractOdeSystem->mNumberOfEquations;
    
    
    std::vector<double> k1(num_equations);
    std::vector<double> k2(num_equations);
    std::vector<double> k3(num_equations);
    std::vector<double> k4(num_equations);
    std::vector<double> yk2(num_equations);
    std::vector<double> yk3(num_equations);
    std::vector<double> yk4(num_equations);
    
    std::vector<double> dy(num_equations);
    std::vector<double> next_y_value(num_equations);
    
    dy = pAbstractOdeSystem->EvaluateYDerivatives(time,currentYValue);
    
    for(int i=0;i<num_equations; i++) 
	{
		k1[i]=timeStep*dy[i];
		yk2[i] = currentYValue[i] + 0.5*k1[i];		
	}
    dy = pAbstractOdeSystem->EvaluateYDerivatives(time+0.5*timeStep,yk2);
	
	for(int i=0;i<num_equations; i++) 
	{
		k2[i]=timeStep*dy[i];
		yk3[i] = currentYValue[i] + 0.5*k2[i];		
	}
    dy = pAbstractOdeSystem->EvaluateYDerivatives(time+0.5*timeStep,yk3);        

	for(int i=0;i<num_equations; i++) 
	{
		k3[i]=timeStep*dy[i];
		yk4[i] = currentYValue[i] + k3[i];		
	}
    dy = pAbstractOdeSystem->EvaluateYDerivatives(time+timeStep,yk4);                
	
	for(int i=0;i<num_equations; i++) 
	{
		k4[i]=timeStep*dy[i];
		next_y_value[i] = currentYValue[i] + (k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;		
	}
		
					
	return next_y_value;
}
