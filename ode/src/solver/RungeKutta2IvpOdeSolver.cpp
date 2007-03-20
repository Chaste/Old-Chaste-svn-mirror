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

void RungeKutta2IvpOdeSolver::CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                                                  double timeStep,
                                                  double time,
                                                  std::vector<double>& currentYValues,
                                                  std::vector<double>& nextYValues)
{
    const unsigned num_equations = pAbstractOdeSystem->GetNumberOfStateVariables();
    
    // Apply Runge-Kutta 2nd Order method for each timestep in AbstractOneStepIvpSolver.
    // Calculates a vector containing the next Y value from the current one for each
    // equation in the system.
    
    std::vector<double> k1(num_equations);
    std::vector<double>& dy = nextYValues; // re-use memory
    
    // Work out k1
    pAbstractOdeSystem->EvaluateYDerivatives(time, currentYValues, dy);
    
    for (unsigned i=0; i<num_equations; i++)
    {
        k1[i] = timeStep*dy[i];
        k1[i] = k1[i]/2.0+currentYValues[i];
    }
    
    // Work out k2 and new solution
    pAbstractOdeSystem->EvaluateYDerivatives(time+timeStep/2.0, k1, dy);
    for (unsigned i=0; i<num_equations; i++)
    {
        nextYValues[i] = currentYValues[i] + timeStep*dy[i];
    }
}
