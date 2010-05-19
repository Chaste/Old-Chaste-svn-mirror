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

#include "OdeSolution.hpp"
#include "PetscTools.hpp"
#include "AbstractOdeSystem.hpp"
#include "Exception.hpp"

OdeSolution::OdeSolution()
   : mNumberOfTimeSteps(0u),
     mpOdeSystemInformation()
{
}

unsigned OdeSolution::GetNumberOfTimeSteps()
{
    return mNumberOfTimeSteps;
}

void OdeSolution::SetNumberOfTimeSteps(unsigned numTimeSteps)
{
    mNumberOfTimeSteps = numTimeSteps;
    mTimes.reserve(numTimeSteps+1);
    mSolutions.reserve(numTimeSteps);
}

void OdeSolution::SetOdeSystemInformation(boost::shared_ptr<const AbstractOdeSystemInformation> pOdeSystemInfo)
{
    mpOdeSystemInformation = pOdeSystemInfo;
}

std::vector<double> OdeSolution::GetVariableAtIndex(unsigned index)
{
    std::vector<double> answer;
    answer.reserve(mSolutions.size());
    for (unsigned i=0; i<mSolutions.size(); i++)
    {
        answer.push_back(mSolutions[i][index]);
    }
    return answer;
}



std::vector<double>& OdeSolution::rGetTimes()
{
    return mTimes;
}

std::vector<std::vector<double> >& OdeSolution::rGetSolutions()
{
    return mSolutions;
}

std::vector<std::vector<double> >& OdeSolution::rGetDerivedQuantities(AbstractOdeSystem* pOdeSystem)
{
    assert(pOdeSystem != NULL);
    if (mDerivedQuantities.empty() && pOdeSystem->GetNumberOfDerivedQuantities() > 0)
    {
        assert(mTimes.size() == mSolutions.size()); // Paranoia
        mDerivedQuantities.reserve(mTimes.size());
        for (unsigned i=0; i<mTimes.size(); i++)
        {
            mDerivedQuantities.push_back(pOdeSystem->ComputeDerivedQuantities(mTimes[i], mSolutions[i]));
        }
    }
    return mDerivedQuantities;
}

void OdeSolution::WriteToFile(std::string directoryName,
                              std::string baseResultsFilename,
                              std::string timeUnits,
                              unsigned stepsPerRow,
                              bool cleanDirectory,
                              unsigned precision,
                              bool includeDerivedQuantities,
                              AbstractOdeSystem* pOdeSystem)
{
    assert(stepsPerRow > 0);
    assert(mTimes.size() > 0);
    assert(mTimes.size() == mSolutions.size());
    assert(mpOdeSystemInformation.get() != NULL);
    if (includeDerivedQuantities)
    {
        if (pOdeSystem == NULL)
        {
            EXCEPTION("You must provide an ODE system to compute derived quantities.");
        }
        assert(pOdeSystem->GetSystemInformation() == mpOdeSystemInformation); // Just in case...
    }

    // Write data to a file using ColumnDataWriter
    ColumnDataWriter writer(directoryName, baseResultsFilename, cleanDirectory, precision);

    if (!PetscTools::AmMaster())
    {
        //Only the master actually writes to file
        return;
    }

    int time_var_id = writer.DefineUnlimitedDimension("Time", timeUnits);

    // Either: the ODE system should have no names&units defined, or it should
    // the same number as the number of solutions per timestep.
    assert( mpOdeSystemInformation->rGetStateVariableNames().size()==0 ||
            (mpOdeSystemInformation->rGetStateVariableNames().size()==mSolutions[0].size()) );

    unsigned num_vars = mSolutions[0].size();

    std::vector<int> var_ids;
    var_ids.reserve(num_vars);
    if (mpOdeSystemInformation->rGetStateVariableNames().size() > 0)
    {
        for (unsigned i=0; i<num_vars; i++)
        {
            var_ids.push_back(writer.DefineVariable(mpOdeSystemInformation->rGetStateVariableNames()[i],
                                                    mpOdeSystemInformation->rGetStateVariableUnits()[i]));
        }
    }
    else
    {
        for (unsigned i=0; i<num_vars; i++)
        {
            std::stringstream string_stream;
            string_stream << "var_" << i;
            var_ids.push_back(writer.DefineVariable(string_stream.str(), ""));
        }
    }
    
    if (includeDerivedQuantities)
    {
        var_ids.reserve(num_vars + pOdeSystem->GetNumberOfDerivedQuantities());
        for (unsigned i=0; i<pOdeSystem->GetNumberOfDerivedQuantities(); i++)
        {
            var_ids.push_back(writer.DefineVariable(mpOdeSystemInformation->rGetDerivedQuantityNames()[i],
                                                    mpOdeSystemInformation->rGetDerivedQuantityUnits()[i]));
        }
    }

    writer.EndDefineMode();

    for (unsigned i=0; i<mSolutions.size(); i+=stepsPerRow)
    {
        writer.PutVariable(time_var_id, mTimes[i]);
        for (unsigned j=0; j<num_vars; j++)
        {
            writer.PutVariable(var_ids[j], mSolutions[i][j]);
        }
        if (includeDerivedQuantities)
        {
            for (unsigned j=0; j<pOdeSystem->GetNumberOfDerivedQuantities(); j++)
            {
                writer.PutVariable(var_ids[j+num_vars], rGetDerivedQuantities(pOdeSystem)[i][j]);
            }
        }
        writer.AdvanceAlongUnlimitedDimension();
    }
    writer.Close();
}
