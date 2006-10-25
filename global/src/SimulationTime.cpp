/**
 * Simulation time object stores the simulation time, uses the singleton pattern
 */

#include "SimulationTime.hpp"
#include "Exception.hpp"
#include <assert.h>

SimulationTime* SimulationTime::mpInstance = 0;

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


void SimulationTime::Destroy()
{
    delete mpInstance;
    mpInstance = 0;
}

/**
 * Get the simlation time step.
 * Warning: Use of this method may result in round errors
 *  -- see GetDimensionalisedTime.
 * @return time step
 */
double SimulationTime::GetTimeStep()
{
    assert(mEndTimeAndNumberOfTimeStepsSet==true);
    return mDurationOfSimulation/mTotalTimeStepsInSimulation;
}

/**
 * Increment the simulation time by one time step
 */
void SimulationTime::IncrementTimeOneStep()
{
    assert(mEndTimeAndNumberOfTimeStepsSet==true);
    mTimeStepsElapsed++;
}

/**
 * Get the number of time steps that have elapsed.
 * @return number of time steps
 */
int SimulationTime::GetTimeStepsElapsed()
{
    assert(mEndTimeAndNumberOfTimeStepsSet==true);
    return mTimeStepsElapsed;
}

/**
 * Get the dimensionalised simulation time.
 * Should not have rounding errors.
 * @return simulation time
 */
double SimulationTime::GetDimensionalisedTime()
{
    assert(mEndTimeAndNumberOfTimeStepsSet==true);
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
void SimulationTime::SetEndTimeAndNumberOfTimeSteps(double durationOfSimulation, int totalTimeStepsInSimulation)
{    
    assert(mEndTimeAndNumberOfTimeStepsSet==false);
    mDurationOfSimulation = durationOfSimulation;
    mTotalTimeStepsInSimulation=totalTimeStepsInSimulation;
    mEndTimeAndNumberOfTimeStepsSet = true;
}
