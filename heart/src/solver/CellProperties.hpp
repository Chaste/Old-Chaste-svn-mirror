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


#ifndef _CELLPROPERTIES_HPP_
#define _CELLPROPERTIES_HPP_

#include "OdeSolution.hpp"

/**
 * Class to calculate various physiological properties from the results of a
 * cardiac simulation.
 *
 * It will calculate for a single cell:
 *   Action potential duration (at any percentage)
 *   Max. upstroke velocity
 *   Action potential amplitude
 *   Max. & min. membrane potential
 *   Cycle length (time between APs)
 * These will all be calculated for the last complete action potential in a
 * dataset.
 *
 * Results may be incorrect if you are stimulating the cell strangely.
 */

class CellProperties
{
private:
    /**
     * The simulation results to process
     */
    std::vector<double>& mrVoltage;
    std::vector<double>& mrTime;

    /**
     * Threshold for determining what counts as an action potential.
     * This is a value part way between the min & max potential, to avoid
     * problems due to 'notches' in an action potential.
     */
    double mThreshold;

    /**
     * Cached vectors containing AP properties
     */
    std::vector<double> mOnsets;
    std::vector<double> mRestingValues;
    std::vector<double> mCycleLengths;
    std::vector<double> mPeakValues;
    std::vector<double> mMaxUpstrokeVelocities;
    std::vector<double> mTimesAtMaxUpstrokeVelocity;

    /**
     * Calculate all the cacheable values.
     */
    void CalculateProperties();

    /**
     * Actually calculate APD.
     * 
     * APD is taken to be the time from the max-upstroke-velocity to 'percentage' of the 
     * 'AP global max potential'--'Potential from just before the start of the AP' range.
     *
     * @param percentage  The percentage of the amplitude to calculate for.
     * @param onset  The time at which the upstroke reaches the threshold.
     * @param minPotential  The minimum potential of this AP.
     * @param maxPotential  The maximum potential of this AP.
     */
    std::vector<double> CalculateActionPotentialDurations(const double percentage,
                                            std::vector<double>& rOnsets,
                                            std::vector<double>& rRestingPotentials,
                                            std::vector<double>& rPeakPotentials);

public:
    /**
     * Constructor sets the data and calls CalculateProperties
     */
    CellProperties(std::vector<double> &rVoltage, std::vector<double> &rTime,  double threshold=-30.0)
        : mrVoltage(rVoltage),
          mrTime(rTime),
          mThreshold(threshold)
    {
        CalculateProperties();
    }


    /**
     * Return the maximum upstroke velocity for all APs.
     */
    std::vector<double> GetMaxUpstrokeVelocities()
    {
        return mMaxUpstrokeVelocities;
    }
    
     /**
     * Return the maximum upstroke velocity for the last AP.
     */
    double GetLastMaxUpstrokeVelocity();
    
    /**
    * Return the time at which the maximum upstroke velocity occured for all APs.
    */
    std::vector<double> GetTimesAtMaxUpstrokeVelocity()
    {
        return mTimesAtMaxUpstrokeVelocity;
    }
    /**
    * Return the time at which the maximum upstroke velocity for the last AP occurred.
    * Returns -1 if no complete AP occurred
    */
    double GetTimeAtLastMaxUpstrokeVelocity();
    
    /**
     * Return the cycle lengths for all APs.
     */
    std::vector<double>  GetCycleLengths()
    {
        return mCycleLengths;
    }
    
    /**
     * Return the peak potentials for all APs.
     */
    std::vector<double>  GetPeakPotentials()
    {
        return mPeakValues;
    }
    /**
     * Return the resting potentials before each AP.
     */
    std::vector<double> GetRestingPotentials()
    {
        return mRestingValues;
    }

    /**
     * Returns all the action potentials
     * 
     * @param percentage is the repolarisation percentage that 
     * the APD will be calculated for. e.g. percentage = 90 for APD90.
     */
    std::vector<double> GetAllActionPotentialDurations(const double percentage);
  
     /**
     * Return the amplitude of the last action potential generated.
     * Throws an exception if no AP is generated.
     * 
     * @param percentage is the repolarisation percentage that 
     * the APD will be calculated for. e.g. percentage = 90 for APD90.
     */  
    double GetLastActionPotentialDuration(const double percentage);
    
     /**
     * Return the amplitude of all the action potentials calculated.
     */
    std::vector<double> GetActionPotentialAmplitudes();
};

#endif //_CELLPROPERTIES_HPP_
