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
#include "SimulationTime.hpp"
#include "Exception.hpp"
#include <cassert>


/** Pointer to the single instance */
SimulationTime* SimulationTime::mpInstance = NULL;

SimulationTime* SimulationTime::Instance()
{
    if (mpInstance == NULL)
    {
        mpInstance = new SimulationTime;
    }
    return mpInstance;
}


SimulationTime::SimulationTime()
    : mTimeStepsElapsed(0),
      mEndTimeAndNumberOfTimeStepsSet(false),
      mCurrentTime(0.0),
      mStartTime(0.0),
      mStartTimeSet(false)
{
    // Make sure there's only one instance - enforces correct serialization
    assert(mpInstance == NULL);
}


void SimulationTime::Destroy()
{
    if (mpInstance)
    {
        delete mpInstance;
        mpInstance = NULL;
    }
}


void SimulationTime::SetStartTime(double startTime)
{
    assert(mStartTimeSet==false);
    mStartTime = startTime;
    mCurrentTime = startTime;
    mStartTimeSet = true;
}


double SimulationTime::GetTimeStep() const
{
    assert(mStartTimeSet);
    assert(mEndTimeAndNumberOfTimeStepsSet);
    return mDurationOfSimulation/mTotalTimeStepsInSimulation;
}


void SimulationTime::IncrementTimeOneStep()
{
    assert(mStartTimeSet);
    assert(mEndTimeAndNumberOfTimeStepsSet);
    mTimeStepsElapsed++;
    mCurrentTime = mStartTime
                   + ((double)mTimeStepsElapsed / (double)mTotalTimeStepsInSimulation)
                   * mDurationOfSimulation;
}


unsigned SimulationTime::GetTimeStepsElapsed() const
{
    assert(mEndTimeAndNumberOfTimeStepsSet);
    return mTimeStepsElapsed;
}


double SimulationTime::GetTime() const
{
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
    // IMPORTANT NOTE: if this assertion fails, it may be because Destroy   //
    // wasn't called in the previous test                                   //
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
    assert(mStartTimeSet);

    return mCurrentTime;
}


void SimulationTime::SetEndTimeAndNumberOfTimeSteps(double endTime, unsigned totalTimeStepsInSimulation)
{
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
    // IMPORTANT NOTE: if this assertion fails, it may be because Destroy   //
    // wasn't called in the previous test                                   //
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
    assert(mStartTimeSet);
    assert(!mEndTimeAndNumberOfTimeStepsSet);
    assert(endTime>mCurrentTime);
    mEndTime = endTime;
    mDurationOfSimulation = mEndTime - mCurrentTime;
    mTotalTimeStepsInSimulation = totalTimeStepsInSimulation;
    mEndTimeAndNumberOfTimeStepsSet = true;
}


void SimulationTime::ResetEndTimeAndNumberOfTimeSteps(const double& rEndTime, const unsigned& rNumberOfTimeStepsInThisRun)
{
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
    // IMPORTANT NOTE: if this assertion fails, it may be because Destroy   //
    // wasn't called in the previous test                                   //
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
    assert(mStartTimeSet);
    assert(rEndTime>mCurrentTime);
    assert(mTimeStepsElapsed>0); // if this throws you should be using set rather than reset.
    //reset the machinery that works out the time
    mStartTime = mCurrentTime;
    mTimeStepsElapsed = 0;
    // set up the new end time and stuff
    mEndTime = rEndTime;
    mDurationOfSimulation = mEndTime - mCurrentTime;
    mTotalTimeStepsInSimulation = rNumberOfTimeStepsInThisRun;
    mEndTimeAndNumberOfTimeStepsSet = true;
}


bool SimulationTime::IsStartTimeSetUp() const
{
    return mStartTimeSet;
}


bool SimulationTime::IsFinished() const
{
    return(mCurrentTime>=mEndTime);
}


unsigned SimulationTime::GetTotalNumberOfTimeSteps() const
{
    return mTotalTimeStepsInSimulation;
}
