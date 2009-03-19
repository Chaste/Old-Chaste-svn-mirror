/*

Copyright (C) University of Oxford, 2005-2009

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
 * OdeSolution.  Sets us the class of ODE solutions, including a function that
 * allows us to save the output data to file.
 */
#ifndef _ODESOLUTION_HPP_
#define _ODESOLUTION_HPP_

#include <vector>
#include <fstream>
#include <cassert>
#include "ColumnDataWriter.hpp"
#include "AbstractOdeSystem.hpp"
#include <sstream>

class OdeSolution
{
private:

    /** Variable for the number of timesteps. */
    unsigned mNumberOfTimeSteps;

    /** A vector of times at each timestep. */
    std::vector<double> mTimes;

    /** Solutions for each variable at each timestep. */
    std::vector<std::vector<double> > mSolutions;

public:

    /**
     * Get the number of timesteps.
     * 
     * @return mNumberOfTimeSteps
     */
    unsigned GetNumberOfTimeSteps(void)
    {
        return mNumberOfTimeSteps;
    }

    /**
     * Get the number of timesteps.
     * 
     * @param numTimeSteps the number of timesteps to use
     */
    void SetNumberOfTimeSteps(unsigned numTimeSteps)
    {
        mNumberOfTimeSteps = numTimeSteps;
        mTimes.reserve(numTimeSteps+1);
        mSolutions.reserve(numTimeSteps);
    }

    /**
     * Get the values of a state variable with a given index in 
     * the ODE system at each timestep.
     * 
     * @param index  the index of the state variable in the system
     */
    std::vector<double> GetVariableAtIndex(unsigned index)
    {
        std::vector<double> answer;
        answer.reserve(mSolutions.size());
        for (unsigned i=0; i<mSolutions.size(); i++)
        {
            answer.push_back(mSolutions[i][index]);
        }
        return answer;
    }

    /**
     * Get the times at which the solution to the ODE system is stored.
     * 
     * @return mTimes.
     */
    std::vector<double>& rGetTimes()
    {
        return mTimes;
    }

    /**
     * Get the values of the solution to the ODE system at each timestep.
     * 
     * @return mSolutions.
     */
    std::vector<std::vector<double> >& rGetSolutions()
    {
        return mSolutions;
    }

    /**
     * Write the data to a file.
     * 
     * @param directoryName  the directory in which to write the data to file
     * @param baseResultsFilename  the name of the file in which to write the data
     * @param pOdeSystem  pointer to the ODE system solved to obtain these results 
     *                     (needed for state variable names and units)
     * @param timeUnites  name of the units of time used
     * @param stepsPerRow  the solution to the ODE system is written to file every 
     *                    this number of timesteps (defaults to 1)
     * @param cleanDirectory  whether to clean the directory (defaults to true)
     */
    void WriteToFile(std::string directoryName,
                     std::string baseResultsFilename,
                     AbstractOdeSystem* pOdeSystem,
                     std::string timeUnits,
                     unsigned stepsPerRow = 1,
                     bool cleanDirectory = true)
    {
        assert(stepsPerRow > 0);
        assert(mTimes.size() > 0);
        assert(mTimes.size() == mSolutions.size());

        // Write data to a file using ColumnDataWriter
        ColumnDataWriter writer(directoryName, baseResultsFilename, cleanDirectory);

        int time_var_id = writer.DefineUnlimitedDimension("Time", timeUnits);

        // Either: the ODE system should have no names&units defined, or it should
        // the same number as the number of solutions per timestep.
        assert( pOdeSystem->rGetVariableNames().size()==0 ||
                (pOdeSystem->rGetVariableNames().size()==mSolutions[0].size()) );

        unsigned num_vars = mSolutions[0].size();

        std::vector<int> var_ids;
        var_ids.reserve(num_vars);
        if (pOdeSystem->rGetVariableNames().size() > 0)
        {
            for (unsigned i=0; i<num_vars; i++)
            {
                var_ids.push_back(writer.DefineVariable(pOdeSystem->rGetVariableNames()[i],
                                                        pOdeSystem->rGetVariableUnits()[i]));
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

        writer.EndDefineMode();

        for (unsigned i=0; i<mSolutions.size(); i+=stepsPerRow)
        {
            writer.PutVariable(time_var_id, mTimes[i]);
            for (unsigned j=0; j<var_ids.size(); j++)
            {
                writer.PutVariable(var_ids[j], mSolutions[i][j]);
            }
            writer.AdvanceAlongUnlimitedDimension();
        }
        writer.Close();
    }
};

#endif //_ODESOLUTION_HPP_
