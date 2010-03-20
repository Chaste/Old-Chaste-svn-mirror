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

#include "AbstractOdeSystem.hpp"

AbstractOdeSystem::AbstractOdeSystem(unsigned numberOfStateVariables)
    : mNumberOfStateVariables(numberOfStateVariables),
      mUseAnalyticJacobian(false)
{
}

AbstractOdeSystem::~AbstractOdeSystem()
{}

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
    const std::vector<std::string>& r_var_names = rGetVariableNames();
    assert(Y.size() == r_var_names.size());
    for (unsigned i=0; i<r_var_names.size(); i++)
    {
        res << "\t" << r_var_names[i] << ":" << Y[i] << "\n";
    }
    return res.str();
}

unsigned AbstractOdeSystem::GetNumberOfStateVariables() const
{
    return mNumberOfStateVariables;
}


unsigned AbstractOdeSystem::GetNumberOfParameters() const
{
    return mParameters.size();
}

double AbstractOdeSystem::GetParameter(unsigned index) const
{
    return mParameters[index];
}

void AbstractOdeSystem::SetParameter(unsigned index, double value)
{
    mParameters[index] = value;
}

const std::vector<std::string>& AbstractOdeSystem::rGetParameterNames() const
{
    assert(mpSystemInfo);
    return mpSystemInfo->rGetParameterNames();
}

const std::vector<std::string>& AbstractOdeSystem::rGetParameterUnits() const
{
    assert(mpSystemInfo);
    return mpSystemInfo->rGetParameterUnits();
}


void AbstractOdeSystem::SetInitialConditions(const std::vector<double>& rInitialConditions)
{
    if (rInitialConditions.size() != mNumberOfStateVariables)
    {
        EXCEPTION("The number of initial conditions must be that of the number of state variables");
    }
    assert(mpSystemInfo);
    mpSystemInfo->SetInitialConditions(rInitialConditions);
}

void AbstractOdeSystem::SetInitialConditionsComponent(unsigned index, double initialCondition)
{
    if (index >= mNumberOfStateVariables)
    {
        EXCEPTION("Index is greater than the number of state variables");
    }
    assert(mpSystemInfo);
    mpSystemInfo->SetInitialConditionsComponent(index, initialCondition);
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
        EXCEPTION("The size of the passed in vector must be that of the number of state variables");
    }
    mStateVariables = rStateVariables;
}

void AbstractOdeSystem::SetStateVariable(unsigned stateVariable, double newValue)
{
    if ( mNumberOfStateVariables <= stateVariable )
    {
        EXCEPTION("The index passed in must be less than the number of state variables");
    }
    mStateVariables[stateVariable] = newValue;
}


std::vector<double>& AbstractOdeSystem::rGetStateVariables()
{
    return mStateVariables;
}

const std::vector<std::string>& AbstractOdeSystem::rGetVariableNames() const
{
    assert(mpSystemInfo);
    return mpSystemInfo->rGetVariableNames();
}

const std::vector<std::string>& AbstractOdeSystem::rGetVariableUnits() const
{
    assert(mpSystemInfo);
    return mpSystemInfo->rGetVariableUnits();
}

boost::shared_ptr<const AbstractOdeSystemInformation> AbstractOdeSystem::GetSystemInformation() const
{
    assert(mpSystemInfo);
    return mpSystemInfo;
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

unsigned AbstractOdeSystem::GetStateVariableNumberByName(const std::string name) const
{
    assert(mpSystemInfo);
    return mpSystemInfo->GetStateVariableNumberByName(name);
}

double AbstractOdeSystem::GetStateVariableValueByNumber(unsigned varNumber) const
{
    assert(varNumber < mNumberOfStateVariables);
    return mStateVariables[varNumber];
}

std::string AbstractOdeSystem::GetStateVariableUnitsByNumber(unsigned varNumber) const
{
    assert(varNumber < mNumberOfStateVariables);
    assert(mpSystemInfo);
    return mpSystemInfo->GetStateVariableUnitsByNumber(varNumber);
}
