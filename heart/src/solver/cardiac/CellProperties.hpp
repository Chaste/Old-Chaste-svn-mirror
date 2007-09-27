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
    std::vector<double> mTime;
    std::vector<double> mVoltage;
    
    /**
     * Threshold for determining what counts as an action potential.
     * This is a value part way between the min & max potential, to avoid
     * problems due to 'notches' in an action potential.
     */
    double mThreshold;
    
    /**
     * Cached values of the properties
     */
    double mMaxUpstrokeVelocity;
    double mTimeAtMaxUpstrokeVelocity;
    double mCycleLength;
    double mMaxPotential, mMinPotential;
    double mUpstrokeStartTime;
    
    /**
     * Values needed for calculating APD.
     */
    double mOnset, mPrevOnset, mPrevMinPotential, mPrevMaxPotential;
    
    /**
     * Calculate all the cacheable values.
     */
    void CalculateProperties();
    
    /**
     * Actually calculate APD.
     * 
     * @param percentage  The percentage of the amplitude to calculate for.
     * @param onset  The time at which the upstroke reaches the threshold.
     * @param minPotential  The minimum potential of this AP.
     * @param maxPotential  The maximum potential of this AP.
     */
    double CalculateActionPotentialDuration(const double percentage,
                                            const double onset,
                                            const double minPotential,
                                            const double maxPotential);
                                            
public:
    /**
     * Constructor does nothing much
     */
    CellProperties(std::vector<double> &rVoltage, std::vector<double> &rTime,  double threshold=-30.0)
    {
        SetData(rVoltage, rTime, threshold);
    }
    
    /**
     * Set the simulation results to process.
     * 
     * @param rSolutionData  The solutions
     * @param vIndex  The index of the membrane potential within the state
     *                variables
     */
    void SetData(std::vector<double> &rVoltage, std::vector<double> &rTime,  double threshold=-30.0)
    {
        mVoltage = rVoltage;
        mTime = rTime;
        mThreshold = threshold;
        CalculateProperties();
    }
    
    /**
     * Return the maximum upstroke velocity.
     */
    double GetMaxUpstrokeVelocity()
    {
        return mMaxUpstrokeVelocity;
    }
    /**
    * Return the time at which the maximum upstroke velocity occured.
    */
    double GetTimeAtMaxUpstrokeVelocity()
    {
        return mTimeAtMaxUpstrokeVelocity;
    }
    /**
     * Return the cycle length.
     */
    double GetCycleLength()
    {
        return mCycleLength;
    }
    /**
     * Return the maximum potential.
     */
    double GetMaxPotential()
    {
        return mMaxPotential;
    }
    /**
     * Return the minimum potential.
     */
    double GetMinPotential()
    {
        return mMinPotential;
    }
    /**
     * Return the amplitude of the action potential.
     */
    double GetActionPotentialAmplitude()
    {
        return mMaxPotential - mMinPotential;
    }
    /**
     * Return the duration of the action potential at a given percentage.
     * 
     * @param percentage  The percentage of the amplitude to calculate for.
     */
    double GetActionPotentialDuration(const double percentage);
};

#endif //_CELLPROPERTIES_HPP_
