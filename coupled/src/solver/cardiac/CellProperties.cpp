#include "CellProperties.hpp"
#include "Exception.hpp"
#include <cmath>

#include <iostream>

enum APPhases { undefined, upstroke, repolarisation };

void CellProperties::CalculateProperties()
{
    // Check we have some suitable data to process
    if (mTime.size() < 1)
    {
        throw Exception("Insufficient time steps to calculate physiological properties.");
    }
    
    // Reset cached properties
    mMaxUpstrokeVelocity = 0.0;
    mCycleLength = 0.0;
    mMaxPotential = 0.0;
    mMinPotential = 0.0;
    mUpstrokeStartTime = 0.0;
    mTimeAtMaxUpstrokeVelocity = -1.0;
    
    double prev_v = mVoltage[0];
    double prev_t = mTime[0];
    double prev_upstroke_vel = 0.0;
    double max_upstroke_vel = 0.0;
    double time_at_max_upstroke_vel = -1.0;
    mOnset = -1.0, mPrevOnset = -1.0;
    APPhases ap_phase = undefined;
    
    int time_steps = mTime.size()-1; //The number of time steps is the number of intervals
    
    for (int i=1; i<=time_steps; i++)
    {
        double v = mVoltage[i];
        double t = mTime[i];

        double upstroke_vel;
    
        if(i==1)
        {   // artificially set initial upstroke_vel to be zero otherwise
            // end up dividing by zero
            upstroke_vel = 0;
        }
        else
        {
            upstroke_vel = (v - prev_v) / (t - prev_t);
        }
                
        switch (ap_phase)
        {
            case undefined:
                if (v <= mThreshold &&
                    prev_upstroke_vel <= 0 && upstroke_vel > 0)
                {
                    // Start of AP, so record minimum V
                    mPrevMinPotential = mMinPotential;
                    mMinPotential = prev_v;
                    
                    // Maximum velocity on this upstroke
                    max_upstroke_vel = upstroke_vel;
                    
                    ap_phase = upstroke;
                }
                break;
            case upstroke:
                if (prev_upstroke_vel >= 0 && upstroke_vel < 0)
                {
                    // End of upstroke, so record max V
                    mPrevMaxPotential = mMaxPotential;
                    mMaxPotential = prev_v;
                    
                    // Store maximum upstroke vel from this upstroke
                    mMaxUpstrokeVelocity = max_upstroke_vel;
                   
                    mTimeAtMaxUpstrokeVelocity = time_at_max_upstroke_vel;
                    
                    ap_phase = repolarisation;
                }
                else 
                {
                    // Update max. upstroke vel
                    if (upstroke_vel > max_upstroke_vel)
                    {
                        max_upstroke_vel = upstroke_vel;
                        time_at_max_upstroke_vel = t;
                    }
                    
                    // Work out cycle length by comparing when we pass
                    // the threshold on successive upstrokes.
                    if (prev_v <= mThreshold && v > mThreshold)
                    {
                        mPrevOnset = mOnset;
                        // Linear interpolation between timesteps
                        mOnset = prev_t + 
                                 (t-prev_t)/(v-prev_v)*(mThreshold-prev_v);
                        // Did we have an earlier upstroke?
                        if (mPrevOnset >= 0)
                        {
                            mCycleLength = mOnset - mPrevOnset;
                        }
                    }
                }
                break;
            case repolarisation:
                if (v < mThreshold)
                {
                    ap_phase = undefined;
                }
                break;
               
        }
        
        prev_v = v;
        prev_t = t;
        prev_upstroke_vel = upstroke_vel;
    }
    
   
}


double CellProperties::CalculateActionPotentialDuration(
                                            const double percentage,
                                            const double onset,
                                            const double minPotential,
                                            const double maxPotential)
{
    double apd = 0.0;
    
    double target_v = 0.01*percentage*(maxPotential-minPotential);
    
    APPhases ap_phase = undefined;
    int time_steps = mTime.size();
    for (int i=0; i<time_steps; i++)
    {
        double t = mTime[i];

        if (ap_phase == undefined && t >= onset)
        {
            ap_phase = upstroke;
        }
        else
        {
            double v = mVoltage[i];
            if (ap_phase == upstroke && fabs(v - maxPotential) < 1e-12)
            {
                ap_phase = repolarisation;
            }
            else if (ap_phase == repolarisation && maxPotential-v >= target_v)
            {
                // We've found the appropriate end time
                apd = t - onset;
                break;
            }
        }
    }
    
    return apd;
}


double CellProperties::GetActionPotentialDuration(const double percentage)
{
    double apd = 0.0;
    // Check we had at least one upstroke
    
    if (mOnset >= 0.0)
    {
        apd = CalculateActionPotentialDuration(percentage, mOnset,
                                               mMinPotential,
                                               mMaxPotential);
        
        if (apd == 0.0 && mPrevOnset >= 0.0)
        {
            // The last action potential is not complete, so try using
            // the previous one (if it exists).
            apd = CalculateActionPotentialDuration(percentage, mPrevOnset,
                                                   mPrevMinPotential,
                                                   mPrevMaxPotential);
        }
    }
    return apd;
}
