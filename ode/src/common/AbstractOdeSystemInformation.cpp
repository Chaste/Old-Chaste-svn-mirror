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


#include <cassert>

#include "AbstractOdeSystemInformation.hpp"
#include "Exception.hpp"

AbstractOdeSystemInformation::AbstractOdeSystemInformation()
    : mInitialised(false)
{
}

AbstractOdeSystemInformation::~AbstractOdeSystemInformation()
{
}

void AbstractOdeSystemInformation::SetInitialConditions(const std::vector<double>& rInitialConditions)
{
    assert(mInitialised);
    mInitialConditions = rInitialConditions;
}

void AbstractOdeSystemInformation::SetInitialConditionsComponent(unsigned index, double initialCondition) throw(std::out_of_range)
{
    assert(mInitialised);
    mInitialConditions.at(index) = initialCondition;
}

std::vector<double> AbstractOdeSystemInformation::GetInitialConditions() const
{
    assert(mInitialised);
    return mInitialConditions;
}

const std::vector<std::string>& AbstractOdeSystemInformation::rGetVariableNames() const
{
    assert(mInitialised);
    return mVariableNames;
}

const std::vector<std::string>& AbstractOdeSystemInformation::rGetVariableUnits() const
{
    assert(mInitialised);
    return mVariableUnits;
}

unsigned AbstractOdeSystemInformation::GetStateVariableNumberByName(const std::string& rName) const
{
    assert(mInitialised);
    unsigned var_number = 0u;
    std::vector<std::string>::const_iterator it = mVariableNames.begin();
    for ( ; it != mVariableNames.end() && *it != rName; ++it, ++var_number);
    if (it == mVariableNames.end())
    {
        EXCEPTION("State variable '" + rName + "' does not exist");
    }
    return var_number;
}

std::string AbstractOdeSystemInformation::GetStateVariableUnitsByNumber(unsigned varNumber) const
{
    assert(mInitialised);
    return mVariableUnits.at(varNumber);
}


const std::vector<std::string>& AbstractOdeSystemInformation::rGetParameterNames() const
{
    assert(mInitialised);
    return mParameterNames;
}

const std::vector<std::string>& AbstractOdeSystemInformation::rGetParameterUnits() const
{
    assert(mInitialised);
    return mParameterUnits;
}
