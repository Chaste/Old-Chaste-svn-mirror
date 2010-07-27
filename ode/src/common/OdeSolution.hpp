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


#ifndef _ODESOLUTION_HPP_
#define _ODESOLUTION_HPP_

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <boost/shared_ptr.hpp>

#include "AbstractOdeSystemInformation.hpp"
class AbstractOdeSystem;

/**
 * An OdeSolution class that that allows us to save the output data to file.
 */
class OdeSolution
{
private:

    /** Variable for the number of timesteps. */
    unsigned mNumberOfTimeSteps;

    /** A vector of times at each timestep. */
    std::vector<double> mTimes;

    /** Solutions for each variable at each timestep. */
    std::vector<std::vector<double> > mSolutions;
    
    /** Derived quantities at each timestep. */
    std::vector<std::vector<double> > mDerivedQuantities;

    /**
     * Information about the concrete ODE system class.
     *
     * Used to get names and units into the file output.
     */
    boost::shared_ptr<const AbstractOdeSystemInformation> mpOdeSystemInformation;

public:

    /**
     * Public constructor - ensures data is empty to start with.
     */
    OdeSolution();

    /**
     * Get the number of timesteps.
     *
     * @return mNumberOfTimeSteps
     */
    unsigned GetNumberOfTimeSteps();

    /**
     * Set the number of timesteps.
     *
     * @param numTimeSteps the number of timesteps to use
     */
    void SetNumberOfTimeSteps(unsigned numTimeSteps);

    /**
     * Set the ODE system information
     *
     * @param pOdeSystemInfo  ODE system information (used to get the names and units of variables).
     */
    void SetOdeSystemInformation(boost::shared_ptr<const AbstractOdeSystemInformation> pOdeSystemInfo);

    /**
     * Get the values of a state variable with a given index in
     * the ODE system at each timestep.
     *
     * @param index  the index of the state variable in the system
     */
    std::vector<double> GetVariableAtIndex(unsigned index);

    /**
     * Get the times at which the solution to the ODE system is stored.
     *
     * @return mTimes.
     */
    std::vector<double>& rGetTimes();

    /**
     * Get the values of the solution to the ODE system at each timestep.
     *
     * @return mSolutions.
     */
    std::vector<std::vector<double> >& rGetSolutions();
    
    /**
     * Get the derived quantities for this ODE system at each timestep.
     * 
     * @param pOdeSystem  the ODE system which was solved to generate this solution object
     */
    std::vector<std::vector<double> >& rGetDerivedQuantities(AbstractOdeSystem* pOdeSystem);

    /**
     * Write the data to a file.
     *
     * @param directoryName  the directory in which to write the data to file
     * @param baseResultsFilename  the name of the file in which to write the data
     * @param timeUnits  name of the units of time used
     * @param stepsPerRow  the solution to the ODE system is written to file every
     *                    this number of timesteps (defaults to 1)
     * @param cleanDirectory  whether to clean the directory (defaults to true)
     * @param precision the precision with which to write the data (i.e. exactly
     *    how many digits to display after the decimal point).  Defaults to 8.
     *    Must be between 2 and 20 (inclusive).
     * @param includeDerivedQuantities  whether to include derived quantities in the output
     * @param pOdeSystem  the ODE system which was solved to generate this solution object
     *    (only used if includeDerivedQuantities=true)
     */
    void WriteToFile(std::string directoryName,
                     std::string baseResultsFilename,
                     std::string timeUnits,
                     unsigned stepsPerRow=1,
                     bool cleanDirectory=true,
                     unsigned precision=8,
                     bool includeDerivedQuantities=false,
                     AbstractOdeSystem* pOdeSystem=NULL);
};

#endif //_ODESOLUTION_HPP_
