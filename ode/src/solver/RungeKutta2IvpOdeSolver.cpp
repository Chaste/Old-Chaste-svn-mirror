/*
Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>.
*/

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
