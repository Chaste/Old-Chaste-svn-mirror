/*

Copyright (C) University of Oxford, 2008

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

/**
 * Return a pointer to the simulation time object.
 * The first time this is called the simulation time object is created.
 */
SimulationTime* SimulationTime::Instance()
{
    if (mpInstance == NULL)
    {
        mpInstance = new SimulationTime;
    }
    return mpInstance;
}

/**
 * Default simulation time constructor
 *
 * Sets up time, you must set the start time,
 * end time and number of time steps before using the object.
 */
SimulationTime::SimulationTime()
{
    // Make sure there's only one instance - enforces correct serialization
    assert(mpInstance == NULL);

    mEndTimeAndNumberOfTimeStepsSet = false;
    mTimeStepsElapsed = 0;
    mCurrentTime = 0.0;
    mStartTime = 0.0;
    mStartTimeSet = false;
}


/**
 * Destroy the current SimulationTime instance.  The next call to
 * Instance will create a new instance, on which
 * SetEndTimeAndNumberOfTimeSteps must be called again to reset time.
 *
 * This method *must* be called before program exit, to avoid a memory
 * leak.
 */
void SimulationTime::Destroy()
{
    if (mpInstance)
    {
        delete mpInstance;
        mpInstance = NULL;
    }
}

/**
 * Set the start time of the simulation
 *
 * @param startTime  The time at which the simulation begins (usually 0.0 hours)
 */
void SimulationTime::SetStartTime(double startTime)
{
    assert(mStartTimeSet==false);
    mStartTime = startTime;
    mCurrentTime = startTime;
    mStartTimeSet = true;
}

/**
 * Get the simulation time step, set in earlier calls.
 * Warning: Use of this method may result in round errors
 *  -- generally use GetTime() instead.
 * @return time step for this run of the simulation
 */
double SimulationTime::GetTimeStep() const
{
    assert(mStartTimeSet);
    assert(mEndTimeAndNumberOfTimeStepsSet);
    return mDurationOfSimulation/mTotalTimeStepsInSimulation;
}

/**
 * Increment the simulation time by one time step.
 *
 * GetTime() will return an updated current time after this call.
 */
void SimulationTime::IncrementTimeOneStep()
{
    assert(mStartTimeSet);
    assert(mEndTimeAndNumberOfTimeStepsSet);
    mTimeStepsElapsed++;
    mCurrentTime = mStartTime
                   + ((double)mTimeStepsElapsed / (double)mTotalTimeStepsInSimulation)
                   * mDurationOfSimulation;
}

/**
 * Get the number of time steps that have elapsed.
 *
 * @return number of time steps which have been taken
 */
unsigned SimulationTime::GetTimeStepsElapsed() const
{
    assert(mEndTimeAndNumberOfTimeStepsSet);
    return mTimeStepsElapsed;
}

/**
 * Get the simulation time (in hours), should not have rounding errors.
 *
 * @return simulation time
 */
double SimulationTime::GetTime() const
{
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
    // IMPORTANT NOTE: if this assertion fails, it may be because Destroy   //
    // wasn't called in the previous test                                   //
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
    assert(mStartTimeSet);

    return mCurrentTime;
}

/**
 * Sets the end time and the number of time steps.
 * This must be called after SetStartTime() but before using any other methods.
 *
 * @param endTime  Time at which to end this run of the simulation
 * @param totalTimeStepsInSimulation  the number of time steps into which the above will be broken
 *
 */
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

/**
 * Reset method for the end time and the number of time steps, to run the simulation
 * further after a first initial run.
 *
 * @param rEndTime The new end time for this simulation (now extended)
 * note that the simulation will run from the current time to this new end time
 * NOT from 0 to this end time.
 * @param rNumberOfTimeStepsInThisRun the number of time steps to
 * split the next run into.
 */
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


/**
 * Allows lower classes to check whether the simulation time class has been set up before using it
 *
 * @return whether the start time of the simulation has been set.
 */
bool SimulationTime::IsStartTimeSetUp() const
{
    return mStartTimeSet;
}

/**
 * @return whether the simulation has finished.
 */
bool SimulationTime::IsFinished() const
{
    return(mCurrentTime>=mEndTime);
}

/**
 * @return the total number of time steps to be taken in this run.
 */
unsigned SimulationTime::GetTotalNumberOfTimeSteps() const
{
    return mTotalTimeStepsInSimulation;
}
