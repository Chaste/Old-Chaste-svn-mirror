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
#include "SimulationTime.hpp"
#include "Warnings.hpp" ///\todo #1885

/** Pointer to the single instance */
SimulationTime* SimulationTime::mpInstance = NULL;

/** Pointer to the delegated (shadow) class \todo #1885*/
TimeStepper* SimulationTime::mpTimeStepper = NULL;

SimulationTime* SimulationTime::Instance()
{
    if (mpInstance == NULL)
    {
        mpInstance = new SimulationTime;
        std::atexit(Destroy);
    }
    return mpInstance;
}

SimulationTime::SimulationTime()
    :
      mTimeStepsElapsed(0),
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
    if (mpTimeStepper)
    {
        delete mpTimeStepper;
        mpTimeStepper = NULL;
    }
}

double SimulationTime::GetTimeStep() const
{
    assert(mStartTimeSet);
    assert(mEndTimeAndNumberOfTimeStepsSet);
    assert(mpTimeStepper);

    assert(fabs(mDurationOfSimulation/mTotalTimeStepsInSimulation - mpTimeStepper->GetIdealTimeStep()) <=DBL_EPSILON);
    return mDurationOfSimulation/mTotalTimeStepsInSimulation;
    ///\todo #1885 return mpTimeStepper->GetIdealTimeStep();
}

void SimulationTime::IncrementTimeOneStep()
{
    assert(mpTimeStepper);
    mpTimeStepper->AdvanceOneTimeStep(); //This can now throw if the end time has been reached
    assert(mStartTimeSet);
    assert(mEndTimeAndNumberOfTimeStepsSet);
    mTimeStepsElapsed++;
    double dt = ((double)mTimeStepsElapsed *mDurationOfSimulation);
    mCurrentTime = mStartTime + dt/(double)mTotalTimeStepsInSimulation;
    assert( fabs(mCurrentTime - mpTimeStepper->GetTime())<1e-9);
}

unsigned SimulationTime::GetTimeStepsElapsed() const
{
    assert(mEndTimeAndNumberOfTimeStepsSet);
    return mTimeStepsElapsed;
}

double SimulationTime::GetTime() const
{
    // NOTE: if this assertion fails, it may be because Destroy() wasn't called in the previous test
    assert(mStartTimeSet);
    //Check if the time stepping has started
    if (mpTimeStepper)
    {
        if (mCurrentTime != mpTimeStepper->GetTime())
        {
            if( mCurrentTime > 1.0)
            {
                assert( fabs(mCurrentTime -mpTimeStepper->GetTime()) < 2*DBL_EPSILON*mCurrentTime);
            }
            else
            {
                assert( fabs(mCurrentTime -mpTimeStepper->GetTime()) < 2*DBL_EPSILON);
            }
        }
        ///\todo #1885 Possibly breaks high-level test: return mpTimeStepper->GetTime();
        return mCurrentTime;
    }
    //If time stepping hasn't started then we are still at start time
    assert(mCurrentTime == mStartTime);
    return mStartTime;
}




void SimulationTime::SetStartTime(double startTime)
{
    assert(!mStartTimeSet);
    mStartTime = startTime;
    mCurrentTime = startTime;
    mStartTimeSet = true;
}

void SimulationTime::SetEndTimeAndNumberOfTimeSteps(double endTime, unsigned totalTimeStepsInSimulation)
{
    // NOTE: if this assertion fails, it may be because Destroy() wasn't called in the previous test
    assert(mStartTimeSet);

    assert(!mEndTimeAndNumberOfTimeStepsSet);
    assert(!mpTimeStepper);
    assert(endTime > mCurrentTime);

    mEndTime = endTime;
    mDurationOfSimulation = mEndTime - mCurrentTime;
    mTotalTimeStepsInSimulation = totalTimeStepsInSimulation;
    mEndTimeAndNumberOfTimeStepsSet = true;

    mpTimeStepper = new TimeStepper(mStartTime, endTime, (endTime-mCurrentTime)/totalTimeStepsInSimulation, true);
}

void SimulationTime::ResetEndTimeAndNumberOfTimeSteps(const double& rEndTime, const unsigned& rNumberOfTimeStepsInThisRun)
{
    // NOTE: if this assertion fails, it may be because Destroy() wasn't called in the previous test
    assert(mStartTimeSet);
    assert(rEndTime > mCurrentTime);

    // NOTE: If this assertion fails, you should be using set rather than reset
    assert(mTimeStepsElapsed > 0);
    assert(mpTimeStepper);

    // Reset the machinery that works out the time
    assert(fabs(mCurrentTime-mpTimeStepper->GetTime())<=DBL_EPSILON);
    //mStartTime = mpTimeStepper->GetTime();
    mStartTime = mCurrentTime;
    mTimeStepsElapsed = 0;

    // Set up the new end time and other member variables
    mEndTime = rEndTime;
    mDurationOfSimulation = mEndTime - mCurrentTime;
    mTotalTimeStepsInSimulation = rNumberOfTimeStepsInThisRun;
    mEndTimeAndNumberOfTimeStepsSet = true;

    delete mpTimeStepper;
    mpTimeStepper = new TimeStepper(mStartTime, rEndTime, (rEndTime-mCurrentTime)/rNumberOfTimeStepsInThisRun, true);
}

bool SimulationTime::IsStartTimeSetUp() const
{
    return mStartTimeSet;
}

bool SimulationTime::IsFinished() const
{
    if(mpTimeStepper->IsTimeAtEnd() != (mCurrentTime>=mEndTime))
    {
        assert (mpTimeStepper->IsTimeAtEnd());
        WARNING("TimeStepper runs past end time in this test");
    }
    //return(mpTimeStepper->IsTimeAtEnd());
    return (mCurrentTime>=mEndTime);
}

unsigned SimulationTime::GetTotalNumberOfTimeSteps() const
{
    return mTotalTimeStepsInSimulation;
}
