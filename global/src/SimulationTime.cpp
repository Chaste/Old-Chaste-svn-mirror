
#include "SimulationTime.hpp"
#include "Exception.hpp"
#include <assert.h>

/** Pointer to the single instance */
SimulationTime* SimulationTime::mpInstance = NULL;

/**
 * Return a pointer to the simulation time object.
 * The first time this is called the simulation time object is created. 
 * */
SimulationTime* SimulationTime::Instance()
{
    if(mpInstance == NULL)
    {
        mpInstance = new SimulationTime;
    }
    return mpInstance;
}

SimulationTime::SimulationTime()
{
    mEndTimeAndNumberOfTimeStepsSet = false;
    mTimeStepsElapsed = 0;
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
 * Get the simlation time step.
 * Warning: Use of this method may result in round errors
 *  -- see GetDimensionalisedTime.
 * @return time step
 */
double SimulationTime::GetTimeStep()
{
    assert(mEndTimeAndNumberOfTimeStepsSet);
    return mDurationOfSimulation/mTotalTimeStepsInSimulation;
}

/**
 * Increment the simulation time by one time step
 */
void SimulationTime::IncrementTimeOneStep()
{
    assert(mEndTimeAndNumberOfTimeStepsSet);
    mTimeStepsElapsed++;
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
    assert(mEndTimeAndNumberOfTimeStepsSet);
    return ((double)mTimeStepsElapsed / (double)mTotalTimeStepsInSimulation)
           * mDurationOfSimulation;
}

/**
 * Sets the end time and the number of time steps.
 * This must be called before any other methods. 
 * @param durationOfSimulation Total dimensionalized time of the simulation
 * @param totalTimeStepsInSimulation the number of time steps into which the above will be broken
 * 
 */
void SimulationTime::SetEndTimeAndNumberOfTimeSteps(double durationOfSimulation, unsigned totalTimeStepsInSimulation)
{    
    assert(!mEndTimeAndNumberOfTimeStepsSet);
    mDurationOfSimulation = durationOfSimulation;
    mTotalTimeStepsInSimulation = totalTimeStepsInSimulation;
    mEndTimeAndNumberOfTimeStepsSet = true;
}


/**
 * Allows lower classes to check whether the simulation time class has been set up before using it
 */
bool SimulationTime::IsSimulationTimeSetUp()
{
	return mEndTimeAndNumberOfTimeStepsSet;
}
