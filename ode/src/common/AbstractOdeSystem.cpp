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

#include "AbstractOdeSystem.hpp"

AbstractOdeSystem::AbstractOdeSystem(unsigned numberOfStateVariables)
{
    mNumberOfStateVariables = numberOfStateVariables;
    mUseAnalytic = false;
}

AbstractOdeSystem::~AbstractOdeSystem()
{}

unsigned AbstractOdeSystem::GetStateVariableNumberByName(const std::string name)
{
    unsigned var_number=0;
    while (var_number != mNumberOfStateVariables &&
           mVariableNames[var_number] != name)
    {
        var_number++;
    }
    if (var_number == mNumberOfStateVariables)
    {
        EXCEPTION("State variable does not exist");
    }
    return var_number;
}

std::string AbstractOdeSystem::DumpState(const std::string& message,
                                         std::vector<double> Y)
{
    std::stringstream res;
    res << message << "\nState:\n";
    if (Y.empty())
    {
        Y = rGetStateVariables();
    }
    assert(Y.size() == mVariableNames.size());
    for (unsigned i=0; i<mVariableNames.size(); i++)
    {
        res << "\t" << mVariableNames[i] << ":" << Y[i] << "\n";
    }
    return res.str();
}
