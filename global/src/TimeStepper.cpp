/*

Copyright (C) University of Oxford, 2005-2011

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
#include <cmath>

#include "TimeStepper.hpp"
#include "Exception.hpp"
#include "MathsCustomFunctions.hpp"

TimeStepper::TimeStepper(double startTime, double endTime, double dt, bool enforceConstantTimeStep, std::vector<double> additionalTimes)
    : mStart(startTime),
      mEnd(endTime),
      mDt(dt),
      mTotalTimeStepsTaken(0),
      mAdditionalTimesReached(0),
      mTime(startTime),
      mEpsilon(DBL_EPSILON)
{
    if (startTime > endTime)
    {
        EXCEPTION("The simulation duration must be positive, not " << endTime-startTime);
    }

    // Remove any additionalTimes entries which fall too close to a time when the stepper would stop anyway
    for (unsigned i=0; i<additionalTimes.size(); i++)
    {
        if (i > 0)
        {
            if (additionalTimes[i-1] >= additionalTimes[i])
            {
                EXCEPTION("The additional times vector should be in ascending numerical order; "
                          "entry " << i << " is less than or equal to entry " << i-1 << ".");
            }
        }

        double time_interval = additionalTimes[i] - startTime;

        // When mDt divides this interval (and the interval is positive) then we are going there anyway
        if (!Divides(mDt, time_interval) && (time_interval > DBL_EPSILON))
        {
            mAdditionalTimes.push_back(additionalTimes[i]);
        }
    }

    /*
     * Note that when mEnd is large then the error of subtracting two numbers of
     * that magnitude is about DBL_EPSILON*mEnd (1e-16*mEnd). When mEnd is small
     * then the error should be around DBL_EPSILON.
     */
    if (mEnd > 1.0)
    {
        mEpsilon = DBL_EPSILON*mEnd;
    }

    // If enforceConstantTimeStep check whether the times are such that we won't have a variable dt
    if (enforceConstantTimeStep)
    {
        double expected_end_time = mStart + mDt*EstimateTimeSteps();

        if (fabs( mEnd - expected_end_time ) > mEpsilon)
        {
            EXCEPTION("TimeStepper estimates non-constant timesteps will need to be used: check timestep "
                      "divides (end_time-start_time) (or divides printing timestep). "
                      "[End time=" << mEnd << "; start=" << mStart << "; dt=" << mDt << "; error="
                      << fabs(mEnd-expected_end_time) << "]");
        }
    }

    mNextTime = CalculateNextTime();
}

double TimeStepper::CalculateNextTime()
{
    double next_time = mStart + (mTotalTimeStepsTaken - mAdditionalTimesReached + 1)*mDt;

    // Does the next time bring us very close to the end time?
    // Note that the inequality in this guard matches the inversion of the guard in the enforceConstantTimeStep
    // calculation of the constructor
    if (mEnd - next_time <= mEpsilon)
    {
        next_time = mEnd;
    }

    if (!mAdditionalTimes.empty())
    {
        if (mAdditionalTimesReached < mAdditionalTimes.size())
        {
            // Does this next step take us very close to, or over, an additional time?
            double next_additional_time = mAdditionalTimes[mAdditionalTimesReached];
            double epsilon = next_additional_time > 1 ? next_additional_time*DBL_EPSILON : DBL_EPSILON;
            if (next_additional_time - next_time <= epsilon)
            {
                next_time = next_additional_time;
                mAdditionalTimesReached++;
            }
        }
    }
    return next_time;
}

void TimeStepper::AdvanceOneTimeStep()
{
    mTotalTimeStepsTaken++;
    if (mTotalTimeStepsTaken == 0)
    {
        EXCEPTION("Time step counter has overflowed.");
    }
    mTime = mNextTime;

    mNextTime = CalculateNextTime();
}

double TimeStepper::GetTime() const
{
    return mTime;
}

double TimeStepper::GetNextTime() const
{
    return mNextTime;
}

double TimeStepper::GetNextTimeStep()
{
    double dt = mDt;

    if (mNextTime == mEnd)
    {
        dt = mEnd - mTime;
    }

    // If the next time or the current time is one of the additional times, the timestep will not be mDt
    if (mAdditionalTimesReached > 0)
    {
        if ((mNextTime == mAdditionalTimes[mAdditionalTimesReached-1]) || (mTime == mAdditionalTimes[mAdditionalTimesReached-1]))
        {
            dt = mNextTime - mTime;
            assert(dt > 0);
        }
    }

    return dt;
}

bool TimeStepper::IsTimeAtEnd() const
{
    return (mTime >= mEnd);
}

unsigned TimeStepper::EstimateTimeSteps() const
{
    return (unsigned) floor((mEnd - mStart)/mDt + 0.5) + mAdditionalTimes.size();
}

unsigned TimeStepper::GetTotalTimeStepsTaken() const
{
    return mTotalTimeStepsTaken;
}

void TimeStepper::ResetTimeStep(double dt)
{
    assert(dt > 0);
    /*
     * The error in subtracting two numbers of the same magnitude is about
     * DBL_EPSILON times that magnitude (we use the sum of the two numbers
     * here as a conservative estimate of their maximum). When both mDt and
     * dt are small then the error should be around DBL_EPSILON.
     */
    double scale = DBL_EPSILON*(mDt + dt);
    if (mDt + dt < 1.0)
    {
        scale = DBL_EPSILON;
    }
    if (fabs(mDt-dt) > scale)
    {
        mDt = dt;
        mStart = mTime;
        mTotalTimeStepsTaken = 0;

        mNextTime = CalculateNextTime();
    }
}