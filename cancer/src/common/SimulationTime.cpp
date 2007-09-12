
#include "SimulationTime.hpp"
#include "Exception.hpp"
#include <assert.h>
//#include <iostream>

/** Pointer to the single instance */
SimulationTime* SimulationTime::mpInstance = NULL;

/**
 * Return a pointer to the simulation time object.
 * The first time this is called the simulation time object is created.
 * */
SimulationTime* SimulationTime::Instance()
{
    if (mpInstance == NULL)
    {
        mpInstance = new SimulationTime;
    }
    return mpInstance;
}

SimulationTime::SimulationTime()
{
    // Make sure there's only one instance - enforces correct serialization
    assert(mpInstance == NULL);
    
    mEndTimeAndNumberOfTimeStepsSet = false;
    mTimeStepsElapsed = 0;
    mCurrentDimensionalisedTime = 0.0;
    mTimeAtEndOfLastRun = 0.0;
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

void SimulationTime::SetStartTime(double startTime)
{
    assert(mStartTimeSet==false);
    mCurrentDimensionalisedTime = startTime;
    mStartTimeSet = true;
}

/**
 * Get the simlation time step.
 * Warning: Use of this method may result in round errors
 *  -- see GetDimensionalisedTime.
 * @return time step
 */
double SimulationTime::GetTimeStep()
{
    assert(mStartTimeSet);
    assert(mEndTimeAndNumberOfTimeStepsSet);
    return mDurationOfSimulation/mTotalTimeStepsInSimulation;
}

/**
 * Increment the simulation time by one time step
 */
void SimulationTime::IncrementTimeOneStep()
{
    assert(mStartTimeSet);
    assert(mEndTimeAndNumberOfTimeStepsSet);
    mTimeStepsElapsed++;
    mCurrentDimensionalisedTime = mTimeAtEndOfLastRun
                                  + ((double)mTimeStepsElapsed / (double)mTotalTimeStepsInSimulation)
                                  * mDurationOfSimulation;
}

/**
 * Get the number of time steps that have elapsed.
 * @return number of time steps
 */
unsigned SimulationTime::GetTimeStepsElapsed()
{
    assert(mEndTimeAndNumberOfTimeStepsSet);
    return mTimeStepsElapsed;
}

/**
 * Get the dimensionalised simulation time.
 * Should not have rounding errors.
 * @return simulation time
 */
double SimulationTime::GetDimensionalisedTime()
{
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
    // IMPORTANT NOTE: if this assertion fails, it may be because Destroy   //
    // wasn't called in the previous test                                   //
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
    assert(mStartTimeSet);
    
    //std::cout << "Current Time = " << mCurrentDimensionalisedTime << "\n" << std::flush;
    
    return mCurrentDimensionalisedTime;
}

/**
 * Sets the end time and the number of time steps.
 * This must be called before any other methods.
 * @param durationOfSimulation Total dimensionalized time of the simulation
 * @param totalTimeStepsInSimulation the number of time steps into which the above will be broken
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
    assert(endTime>mCurrentDimensionalisedTime);
    mEndTime = endTime;
    mDurationOfSimulation = mEndTime - mCurrentDimensionalisedTime;
    mTotalTimeStepsInSimulation = totalTimeStepsInSimulation;
    mEndTimeAndNumberOfTimeStepsSet = true;
}

/**
 * Reset method for the end time and the number of time steps.
 *
 * @param rEndTime The new end time for this simulation (probably now extended)
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
    assert(rEndTime>mCurrentDimensionalisedTime);
    assert(mTimeStepsElapsed>0); // if this throws you should be using set rather than reset.
    //reset the machinery that works out the time
    mTimeAtEndOfLastRun = mCurrentDimensionalisedTime;
    mTimeStepsElapsed = 0;
    // set up the new end time and stuff
    mEndTime = rEndTime;
    mDurationOfSimulation = mEndTime - mCurrentDimensionalisedTime;
    mTotalTimeStepsInSimulation = rNumberOfTimeStepsInThisRun;
    mEndTimeAndNumberOfTimeStepsSet = true;
}


/**
 * Allows lower classes to check whether the simulation time class has been set up before using it
 */
bool SimulationTime::IsStartTimeSetUp()
{
    return mStartTimeSet;
}

bool SimulationTime::IsFinished()
{
    return(mCurrentDimensionalisedTime>=mEndTime);
}

unsigned SimulationTime::GetTotalNumberOfTimeSteps()
{
    return mTotalTimeStepsInSimulation;
}
