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

#include "TimeStepper.hpp"
#include "Exception.hpp"
#include <cmath>
#include <cfloat>
#include <cassert>

const double SMIDGE = 2e-10;  ///\todo #1827 Correlate this code with the use of Divides in UblasCustomFunctions (which is used for setting time steps)

TimeStepper::TimeStepper(double startTime, double endTime, double dt, bool enforceConstantTimeStep, std::vector<double> additionalTimes)
    : mStart(startTime),
      mEnd(endTime),
      mDt(dt),
      mTotalTimeStepsTaken(0),
      mAdditionalTimesReached(0),
      mTime(startTime)
{
    if (startTime > endTime)
    {
        EXCEPTION("The simulation duration must be positive");
    }

    // Remove any additionalTimes entries which fall too close to a time when the stepper would stop anyway
    for (unsigned i=0;i<additionalTimes.size();i++)
    {
        if (i > 0)
        {
            if (additionalTimes[i-1] >= additionalTimes[i])
            {
                EXCEPTION("The additional times vector should be in ascending numerical order");
            }
        }

        double test_value = (additionalTimes[i]-startTime)/mDt;
        if (fabs(floor(test_value+0.5)-test_value) > 1e-12)
        {
            mAdditionalTimes.push_back(additionalTimes[i]);
        }
    }

    mNextTime = CalculateNextTime();

    // If enforceConstantTimeStep check whether the times are such that we won't have a variable dt
    if (enforceConstantTimeStep)
    {
        if ( fabs(mDt*EstimateTimeSteps()-mEnd+mStart) > SMIDGE )
        {
            EXCEPTION("TimeStepper estimate non-constant timesteps will need to be used: check timestep divides (end_time-start_time) (or divides printing timestep)");
        }
    }
}

double TimeStepper::CalculateNextTime()
{
    double next_time = mStart + (mTotalTimeStepsTaken-mAdditionalTimesReached+1) * mDt;

    if ((next_time) + SMIDGE*(mDt) >= mEnd)
    {
        next_time = mEnd;
    }

    if (    (mAdditionalTimes.size() > 0)                          // any additional times given
         && (mAdditionalTimesReached < mAdditionalTimes.size())    // not yet done all the additional times
         && ((next_time) + SMIDGE*(mDt) >= mAdditionalTimes[mAdditionalTimesReached]) ) // this next step takes us over an additional time
    {
        next_time = mAdditionalTimes[mAdditionalTimesReached];
        mAdditionalTimesReached++;
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
    if ((mAdditionalTimesReached > 0)  &&
        ( (mNextTime == mAdditionalTimes[mAdditionalTimesReached-1]) || (mTime == mAdditionalTimes[mAdditionalTimesReached-1]) ) )
    {
        dt = mNextTime - mTime;
        assert(dt > 0);
    }

    return dt;
}

bool TimeStepper::IsTimeAtEnd() const
{
    return mTime >= mEnd;
}

unsigned TimeStepper::EstimateTimeSteps() const
{
    return (unsigned) floor((mEnd - mStart)/mDt+0.5) + mAdditionalTimes.size();
}

unsigned TimeStepper::GetTotalTimeStepsTaken() const
{
    return mTotalTimeStepsTaken;
}

void TimeStepper::ResetTimeStep(double dt)
{
    if (fabs(mDt-dt) > 1e-8)
    {
        assert(dt > 0);
        mDt = dt;
        mStart = mTime;
        mTotalTimeStepsTaken = 0;

        mNextTime = CalculateNextTime();
    }
}
