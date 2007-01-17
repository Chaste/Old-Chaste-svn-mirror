/**
 * Concrete EulerIvpOdeSolver class.
 */
#include "EulerIvpOdeSolver.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSolution.hpp"

//#include <iostream>
#include <vector>

/**
 * Solves a system of ODEs using the Forward Euler method
 *
 * To be used in the form:
 *
 * EulerIvpOdeSolver mySolver;
 *
 * OdeSolution solution=mySolver.Solve(pMyOdeSystem, yInit, StartTime, EndTime, TimeStep, SamplingTime);
 *
 * See documentation for AbstractOneStepIvpOdeSolver::Solve()
 *
 */

void EulerIvpOdeSolver::CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                                            double timeStep,
                                            double time,
                                            std::vector<double>& currentYValues,
                                            std::vector<double>& nextYValues)
{
    // for each timestep in AbstractOneStepIvpSolver calculates a vector containing
    // the next Y value from the current one for each equation in the system.
    
    const unsigned num_equations = pAbstractOdeSystem->GetNumberOfStateVariables();
    
    // Yes, this looks wierd, but it makes good use of memory!
    pAbstractOdeSystem->EvaluateYDerivatives(time, currentYValues, nextYValues);
    
    for (unsigned i=0;i<num_equations; i++)
    {
        // nextYValues contains dY/dt until here
        nextYValues[i] = currentYValues[i] + timeStep*nextYValues[i];
    }
}
