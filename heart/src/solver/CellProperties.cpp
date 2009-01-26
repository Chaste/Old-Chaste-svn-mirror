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


#include "CellProperties.hpp"
#include "Exception.hpp"
#include <cmath>

#include <iostream>

enum APPhases { UNDEFINED, UPSTROKE, REPOLARISATION };

void CellProperties::CalculateProperties()
{
    // Check we have some suitable data to process
    if (mrTime.size() < 1)
    {
        EXCEPTION("Insufficient time steps to calculate physiological properties.");
    }

    // Reset cached properties
    mMaxUpstrokeVelocity = 0.0;
    mCycleLength = 0.0;
    mMaxPotential = 0.0;
    mMinPotential = 0.0;
    mUpstrokeStartTime = 0.0;
    mTimeAtMaxUpstrokeVelocity = -1.0;

    double prev_v = mrVoltage[0];
    double prev_t = mrTime[0];
    double prev_upstroke_vel = 0.0;
    double max_upstroke_vel = 0.0;
    double time_at_max_upstroke_vel = -1.0;
    mOnset = -1;

    mPrevOnset = -1.0;
    APPhases ap_phase = UNDEFINED;

    unsigned time_steps = mrTime.size()-1; //The number of time steps is the number of intervals

    double current_max_voltage = -DBL_MAX;
    for (unsigned i=1; i<=time_steps; i++)
    {
        double v = mrVoltage[i];
        double t = mrTime[i];

        // calculate the max voltage for this action potential:
        if(v > current_max_voltage)
        {
            current_max_voltage = v;
        }
        if( v<mThreshold && prev_v >= mThreshold )
        {
            // crossing threshold value, so assume ending action potential, so save values
            mPrevMaxPotential = mMaxPotential;
            mMaxPotential = current_max_voltage;
            current_max_voltage = -DBL_MAX;
        }

        double upstroke_vel;

        if (i==1)
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
            case UNDEFINED:
                // First AP: switch on if below threshold and there is a positive gradient
                // Subsequent APs: switch on if below threshold and there is a gradient of at least 5% of the previously seen activation
                // This avoids switching to UPSTROKE if there's a kink on the way down
                if (v <= mThreshold &&  upstroke_vel > 0.05*mMaxUpstrokeVelocity)
                {   
                    
                    // Start of AP, so record minimum V
                    mPrevMinPotential = mMinPotential;
                    mMinPotential = prev_v;

                    // Maximum velocity on this upstroke
                    max_upstroke_vel = upstroke_vel;

                    ap_phase = UPSTROKE;
                }
                break;

            case UPSTROKE:
                if (prev_upstroke_vel >= 0 && upstroke_vel < 0)
                {
                    // Store maximum upstroke vel from this upstroke
                    mMaxUpstrokeVelocity = max_upstroke_vel;
                    mTimeAtMaxUpstrokeVelocity = time_at_max_upstroke_vel;

                    ap_phase = REPOLARISATION;
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

            case REPOLARISATION:
                if (v < mThreshold)
                {
                    ap_phase = UNDEFINED;
                }

                // avoid kinks in the upstroke: if the velocity starts increasing again,
                // go back to upstroke state.
                if(prev_upstroke_vel <=0 && upstroke_vel > 0)
                {
                    ap_phase = UPSTROKE;
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

    APPhases ap_phase = UNDEFINED;
    unsigned time_steps = mrTime.size();
    for (unsigned i=0; i<time_steps; i++)
    {
        double t = mrTime[i];

        if (ap_phase == UNDEFINED && t >= onset)
        {
            ap_phase = UPSTROKE;
        }
        else
        {
            double v = mrVoltage[i];
            if (ap_phase == UPSTROKE && fabs(v - maxPotential) < 1e-12)
            {
                ap_phase = REPOLARISATION;
            }
            else if (ap_phase == REPOLARISATION && maxPotential-v >= target_v)
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
    if (mOnset < 0)
    {
        #define COVERAGE_IGNORE
        // possible false error here if the simulation started at time < 0
        EXCEPTION("No action potential occured");
        #undef COVERAGE_IGNORE
    }

    double apd = CalculateActionPotentialDuration(percentage, mOnset,
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

    return apd;
}
