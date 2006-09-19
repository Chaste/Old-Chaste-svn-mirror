/**
 * Concrete RungeKutta2IvpOdeSolver class.
 */
#include "RungeKutta2IvpOdeSolver.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"

//#include <iostream>
#include <vector>
#include <cassert>

/**
 * Solves a system of ODEs using the Runge Kutta 2nd Order Initial Value Problem Ordinary Differential Equation Solver
 *
 * To be used in the form:
 *
 * RungeKutta2IvpOdeSolver mySolver
 * 
 * OdeSolution solution=mySolver.Solve(pMyOdeSystem, yInit, StartTime, EndTime, TimeStep, SamplingTime);
 * 
 * See documentation for AbstractOneStepIvpOdeSolver::Solve() 
 */

std::vector<double> RungeKutta2IvpOdeSolver::CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
        double timeStep,
        double time,
        std::vector<double> currentYValue)
{
    int num_equations = pAbstractOdeSystem->GetNumberOfStateVariables();
    
    // Apply Runge-Kutta 2nd Order method for each timestep in AbstractOneStepIvpSolver.
    // Calculates a vector containing the next Y value from the current one for each 
    // equation in the system.
    
    // Work out k1
    std::vector<double> k1(num_equations);
    std::vector<double> dy(num_equations);
    std::vector<double> next_y_value(num_equations);
    dy = pAbstractOdeSystem->EvaluateYDerivatives(time,currentYValue);
    for (int i=0;i<num_equations; i++)
    {
        k1[i] = timeStep*dy[i];
        k1[i] = k1[i]/2.0+currentYValue[i];
    }
    // Work out k2
    std::vector<double> k2(num_equations);
    dy = pAbstractOdeSystem->EvaluateYDerivatives(time+timeStep/2.0,k1);
    for (int i=0;i<num_equations; i++)
    {
        k2[i] = timeStep*dy[i];
        
        next_y_value[i]=currentYValue[i] + k2[i];
    }
    
    return next_y_value;
}
