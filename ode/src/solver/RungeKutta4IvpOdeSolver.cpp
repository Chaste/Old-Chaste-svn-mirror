/**
 * Concrete RungeKutta4IvpOdeSolver class.
 */
#include "RungeKutta4IvpOdeSolver.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"

//#include <iostream>
#include <vector>

/**
 * Solves a system of ODEs using the Runge Kutta 4th Order Initial Value Problem Ordinary Differential Equation Solver
 *
 * To be used in the form:
 *
 * RungeKutta4IvpOdeSolver mySolver
 *
 * OdeSolution solution=mySolver.Solve(pMyOdeSystem, yInit, StartTime, EndTime, TimeStep, SamplingTime);
 *
 * See documentation for AbstractOneStepIvpOdeSolver::Solve()
 */

void RungeKutta4IvpOdeSolver::CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
        double timeStep,
        double time,
        std::vector<double>& currentYValues,
        std::vector<double>& nextYValues)
{
    // Apply Runge-Kutta 4th Order method for each timestep in AbstractOneStepIvpSolver.
    // Calculates a vector containing the next Y value from the current one for each
    // equation in the system.
    
    const unsigned num_equations = pAbstractOdeSystem->GetNumberOfStateVariables();
    
    std::vector<double> k1(num_equations);
    std::vector<double> k2(num_equations);
    std::vector<double> k3(num_equations);
    std::vector<double> k4(num_equations);
    std::vector<double> yk2(num_equations);
    std::vector<double> yk3(num_equations);
    std::vector<double> yk4(num_equations);
    
    std::vector<double>& dy = nextYValues; // re-use memory (not that it makes much difference here!)
    
    pAbstractOdeSystem->EvaluateYDerivatives(time, currentYValues, dy);
    for (unsigned i=0;i<num_equations; i++)
    {
        k1[i]=timeStep*dy[i];
        yk2[i] = currentYValues[i] + 0.5*k1[i];
    }
    
    pAbstractOdeSystem->EvaluateYDerivatives(time+0.5*timeStep, yk2, dy);
    for (unsigned i=0;i<num_equations; i++)
    {
        k2[i]=timeStep*dy[i];
        yk3[i] = currentYValues[i] + 0.5*k2[i];
    }
    
    pAbstractOdeSystem->EvaluateYDerivatives(time+0.5*timeStep, yk3, dy);
    for (unsigned i=0;i<num_equations; i++)
    {
        k3[i]=timeStep*dy[i];
        yk4[i] = currentYValues[i] + k3[i];
    }
    
    pAbstractOdeSystem->EvaluateYDerivatives(time+timeStep, yk4, dy);
    for (unsigned i=0;i<num_equations; i++)
    {
        k4[i]=timeStep*dy[i];
        nextYValues[i] = currentYValues[i] + (k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;
    }
}
