#ifndef _PHYSIOLOGICALPROPERTIES_HPP_
#define _PHYSIOLOGICALPROPERTIES_HPP_

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

class PhysiologicalProperties
{
private:
    /**
     * The simulation results to process
     */
    OdeSolution mSolutionData;
    /**
     * The index of the membrane potential within the state variables
     */
    int mVIndex;

    /**
     * Whether we have already calculated the properties and cached the values,
     * or not.
     */
    bool mCalculatedProperties;

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
    PhysiologicalProperties() : mCalculatedProperties(false)
    {
        mSolutionData = OdeSolution();
        mSolutionData.SetNumberOfTimeSteps(0);
    }
    
    /**
     * Set the simulation results to process.
     * 
     * @param rSolutionData  The solutions
     * @param vIndex  The index of the membrane potential within the state
     *                variables
     */
    void SetData(OdeSolution &rSolutionData, int vIndex, double threshold=-30.0)
    {
        mSolutionData = rSolutionData; // Check if does copy or ref
        mVIndex = vIndex;
        mThreshold = threshold;
        mCalculatedProperties = false;
    }
    
    /**
     * Return the maximum upstroke velocity.
     */
    double GetMaxUpstrokeVelocity()
    {
        CalculateProperties();
        return mMaxUpstrokeVelocity;
    }
    /**
     * Return the cycle length.
     */
    double GetCycleLength()
    {
        CalculateProperties();
        return mCycleLength;
    }
    /**
     * Return the maximum potential.
     */
    double GetMaxPotential()
    {
        CalculateProperties();
        return mMaxPotential;
    }
    /**
     * Return the minimum potential.
     */
    double GetMinPotential()
    {
        CalculateProperties();
        return mMinPotential;
    }
    /**
     * Return the amplitude of the action potential.
     */
    double GetActionPotentialAmplitude()
    {
        CalculateProperties();
        return mMaxPotential - mMinPotential;
    }
    /**
     * Return the duration of the action potential at a given percentage.
     * 
     * @param percentage  The percentage of the amplitude to calculate for.
     */
    double GetActionPotentialDuration(const double percentage);
};

#endif //_PHYSIOLOGICALPROPERTIES_HPP_
