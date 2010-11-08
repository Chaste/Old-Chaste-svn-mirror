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
#include <algorithm>

#include "AbstractOdeSystemInformation.hpp"
#include "Exception.hpp"

AbstractOdeSystemInformation::AbstractOdeSystemInformation()
    : mInitialised(false)
{
}

AbstractOdeSystemInformation::~AbstractOdeSystemInformation()
{
}

std::string AbstractOdeSystemInformation::GetSystemName() const
{
    return mSystemName;
}

void AbstractOdeSystemInformation::SetDefaultInitialConditions(const std::vector<double>& rInitialConditions)
{
    assert(mInitialised);
    mInitialConditions = rInitialConditions;
}

void AbstractOdeSystemInformation::SetDefaultInitialCondition(unsigned index, double initialCondition)
{
    assert(mInitialised);
    mInitialConditions.at(index) = initialCondition;
}

std::vector<double> AbstractOdeSystemInformation::GetInitialConditions() const
{
    assert(mInitialised);
    return mInitialConditions;
}

const std::vector<std::string>& AbstractOdeSystemInformation::rGetStateVariableNames() const
{
    assert(mInitialised);
    return mVariableNames;
}

const std::vector<std::string>& AbstractOdeSystemInformation::rGetStateVariableUnits() const
{
    assert(mInitialised);
    return mVariableUnits;
}

unsigned AbstractOdeSystemInformation::GetStateVariableIndex(const std::string& rName) const
{
    assert(mInitialised);
    std::vector<std::string>::const_iterator it = find(mVariableNames.begin(), mVariableNames.end(), rName);
    if (it == mVariableNames.end())
    {
        EXCEPTION("No state variable named '" + rName + "'.");
    }
    return (unsigned)(it - mVariableNames.begin());
}

bool AbstractOdeSystemInformation::HasStateVariable(const std::string& rName) const
{
    assert(mInitialised);
    std::vector<std::string>::const_iterator it = find(mVariableNames.begin(), mVariableNames.end(), rName);
    return (it != mVariableNames.end());
}

std::string AbstractOdeSystemInformation::GetStateVariableUnits(unsigned index) const
{
    assert(mInitialised);
    if (index >= mVariableUnits.size())
    {
        EXCEPTION("The index passed in must be less than the number of state variables.");
    }
    return mVariableUnits[index];
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

unsigned AbstractOdeSystemInformation::GetParameterIndex(const std::string& rName) const
{
    assert(mInitialised);
    std::vector<std::string>::const_iterator it = find(mParameterNames.begin(), mParameterNames.end(), rName);
    if (it == mParameterNames.end())
    {
        EXCEPTION("No parameter named '" + rName + "'.");
    }
    return (unsigned)(it - mParameterNames.begin());
}

bool AbstractOdeSystemInformation::HasParameter(const std::string& rName) const
{
    assert(mInitialised);
    std::vector<std::string>::const_iterator it = find(mParameterNames.begin(), mParameterNames.end(), rName);
    return (it != mParameterNames.end());
}

std::string AbstractOdeSystemInformation::GetParameterUnits(unsigned index) const
{
    assert(mInitialised);
    if (index >= mParameterUnits.size())
    {
        EXCEPTION("The index passed in must be less than the number of parameters.");
    }
    return mParameterUnits[index];
}

unsigned AbstractOdeSystemInformation::GetAnyVariableIndex(const std::string& rName) const
{
    assert(mInitialised);
    if (HasStateVariable(rName))
    {
        return GetStateVariableIndex(rName);
    }
    else if (HasParameter(rName))
    {
        return mVariableNames.size() + GetParameterIndex(rName);
    }
    else if (HasDerivedQuantity(rName))
    {
        return mVariableNames.size() + mParameterNames.size() + GetDerivedQuantityIndex(rName);
    }
    else
    {
        EXCEPTION("No state variable, parameter, or derived quantity named '" + rName + "'.");
    }
}


bool AbstractOdeSystemInformation::HasAnyVariable(const std::string& rName) const
{
    assert(mInitialised);
    return (HasStateVariable(rName) || HasParameter(rName) || HasDerivedQuantity(rName));
}

std::string AbstractOdeSystemInformation::GetAnyVariableUnits(unsigned index) const
{
    assert(mInitialised);
    if (index < mVariableUnits.size())
    {
        return mVariableUnits[index];
    }
    else
    {
        unsigned offset = mVariableUnits.size();
        if (index - offset < mParameterUnits.size())
        {
            return mParameterUnits[index - offset];
        }
        else
        {
            offset += mParameterUnits.size();
            if (index - offset < mDerivedQuantityUnits.size())
            {
                return mDerivedQuantityUnits[index - offset];
            }
            else
            {
                EXCEPTION("Invalid index passed to GetAnyVariableUnits.");
            }
        }
    }
}


const std::vector<std::string>& AbstractOdeSystemInformation::rGetDerivedQuantityNames() const
{
    assert(mInitialised);
    return mDerivedQuantityNames;
}

const std::vector<std::string>& AbstractOdeSystemInformation::rGetDerivedQuantityUnits() const
{
    assert(mInitialised);
    return mDerivedQuantityUnits;
}

unsigned AbstractOdeSystemInformation::GetDerivedQuantityIndex(const std::string& rName) const
{
    assert(mInitialised);
    std::vector<std::string>::const_iterator it = find(mDerivedQuantityNames.begin(), mDerivedQuantityNames.end(), rName);
    if (it == mDerivedQuantityNames.end())
    {
        EXCEPTION("No derived quantity named '" + rName + "'.");
    }
    return (unsigned)(it - mDerivedQuantityNames.begin());
}

bool AbstractOdeSystemInformation::HasDerivedQuantity(const std::string& rName) const
{
    assert(mInitialised);
    std::vector<std::string>::const_iterator it = find(mDerivedQuantityNames.begin(), mDerivedQuantityNames.end(), rName);
    return (it != mDerivedQuantityNames.end());
}

std::string AbstractOdeSystemInformation::GetDerivedQuantityUnits(unsigned index) const
{
    assert(mInitialised);
    if (index >= mDerivedQuantityUnits.size())
    {
        EXCEPTION("The index passed in must be less than the number of derived quantities.");
    }
    return mDerivedQuantityUnits[index];
}


unsigned AbstractOdeSystemInformation::GetNumberOfAttributes() const
{
    assert(mInitialised);
    return mAttributes.size();
}

bool AbstractOdeSystemInformation::HasAttribute(const std::string& rName) const
{
    assert(mInitialised);
    return (mAttributes.find(rName) != mAttributes.end());
}

double AbstractOdeSystemInformation::GetAttribute(const std::string& rName) const
{
    assert(mInitialised);
    std::map<std::string, double>::const_iterator it = mAttributes.find(rName);
    if (it == mAttributes.end())
    {
        EXCEPTION("No attribute '" + rName + "' found.");
    }
    return it->second;
}
