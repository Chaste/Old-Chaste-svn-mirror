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


#include "TimeStepper.hpp"
#include "Exception.hpp"
#include <cmath>
#include <cfloat>
#include <cassert>

const double SMIDGE = 1e-10;

TimeStepper::TimeStepper(double startTime, double endTime, double dt, bool enforceConstantTimeStep)
    : mStart(startTime),
      mEnd(endTime),
      mDt(dt),
      mTimeStep(0),
      mTime(startTime)
{
    if (startTime > endTime)
    {
        EXCEPTION("The simulation duration must be positive");
    }

    /**
     * \todo This assertion breaks several tests
     *
    if (endTime-startTime < dt-SMIDGE)
    {
       std::cout<<"Span is "<<endTime-startTime<<"\n";
       std::cout<<"Delta is "<<dt<<"\n";
       assert(0);
    }
     */

    mNextTime = CalculateNextTime();

    // if enforceConstantTimeStep check whether the times are such that we won't have a variable dt
    if (enforceConstantTimeStep)
    {
        if ( fabs(mDt*EstimateTimeSteps()-mEnd+mStart) > SMIDGE )
        {
            //PRINT_4_VARIABLES(mDt, EstimateTimeSteps(), mDt*EstimateTimeSteps(), mEnd-mStart);
            EXCEPTION("TimeStepper estimate non-constant timesteps will need to be used: check timestep divides (end_time-start_time) (or divides printing timestep)");
        }
    }
}

double TimeStepper::CalculateNextTime() const
{
    double next_time = mStart + (mTimeStep+1) * mDt;
    if ((next_time) + SMIDGE*(mDt) >= mEnd)
    {
        next_time = mEnd;
    }
    return next_time;
}

void TimeStepper::AdvanceOneTimeStep()
{
    mTimeStep++;
    if (mTimeStep == 0)
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

    return dt;
}

bool TimeStepper::IsTimeAtEnd() const
{
    return mTime >= mEnd;
}

unsigned TimeStepper::EstimateTimeSteps() const
{
    return (unsigned) round((mEnd - mStart)/mDt);
}

unsigned TimeStepper::GetTimeStepsElapsed() const
{
    return mTimeStep;
}
