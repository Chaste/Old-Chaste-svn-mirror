/*

Copyright (C) University of Oxford, 2008

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


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

    // Yes, this looks weird, but it makes good use of memory!
    pAbstractOdeSystem->EvaluateYDerivatives(time, currentYValues, nextYValues);

    for (unsigned i=0;i<num_equations; i++)
    {
        // nextYValues contains dY/dt until here
        nextYValues[i] = currentYValues[i] + timeStep*nextYValues[i];
    }
}
