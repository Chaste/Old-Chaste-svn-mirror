/*

Copyright (C) University of Oxford, 2005-2010

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

#include <sstream>
#include <cassert>

#include "AbstractOdeSystem.hpp"
#include "Exception.hpp"

AbstractOdeSystem::AbstractOdeSystem(unsigned numberOfStateVariables)
    : AbstractParameterisedSystem<std::vector<double> >(numberOfStateVariables),
      mUseAnalyticJacobian(false)
{
}

AbstractOdeSystem::~AbstractOdeSystem()
{
}

bool AbstractOdeSystem::CalculateStoppingEvent(double time, const std::vector<double>& rY)
{
    return false;
}

std::string AbstractOdeSystem::DumpState(const std::string& rMessage,
                                         std::vector<double> Y)
{
    std::stringstream res;
    res << rMessage << "\nState:\n";
    if (Y.empty())
    {
        Y = rGetStateVariables();
    }
    const std::vector<std::string>& r_var_names = rGetStateVariableNames();
    assert(Y.size() == r_var_names.size());
    for (unsigned i=0; i<r_var_names.size(); i++)
    {
        res << "\t" << r_var_names[i] << ":" << Y[i] << "\n";
    }
    return res.str();
}


void AbstractOdeSystem::SetInitialConditions(const std::vector<double>& rInitialConditions)
{
    if (rInitialConditions.size() != mNumberOfStateVariables)
    {
        EXCEPTION("The number of initial conditions must be that of the number of state variables.");
    }
    assert(mpSystemInfo);
    mpSystemInfo->SetInitialConditions(rInitialConditions);
}

void AbstractOdeSystem::SetInitialCondition(unsigned index, double initialCondition)
{
    if (index >= mNumberOfStateVariables)
    {
        EXCEPTION("Index is greater than the number of state variables.");
    }
    assert(mpSystemInfo);
    mpSystemInfo->SetInitialCondition(index, initialCondition);
}

std::vector<double> AbstractOdeSystem::GetInitialConditions() const
{
    assert(mpSystemInfo);
    return mpSystemInfo->GetInitialConditions();
}

void AbstractOdeSystem::SetStateVariables(const std::vector<double>& rStateVariables)
{
    if ( mNumberOfStateVariables != rStateVariables.size() )
    {
        EXCEPTION("The size of the passed in vector must be that of the number of state variables.");
    }
    mStateVariables = rStateVariables;
}

std::vector<double>& AbstractOdeSystem::rGetStateVariables()
{
    return mStateVariables;
}


double AbstractOdeSystem::CalculateRootFunction(double time, const std::vector<double>& rY)
{
    bool stop = CalculateStoppingEvent(time, rY);
    return stop ? 0.0 : 1.0;
}

bool AbstractOdeSystem::GetUseAnalyticJacobian()
{
    return mUseAnalyticJacobian;
}
