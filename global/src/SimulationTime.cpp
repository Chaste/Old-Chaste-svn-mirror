/**
 * Simulation time object stores the simulation time, uses the singleton pattern
 */

#include "SimulationTime.hpp"
#include "Exception.hpp"
#include <assert.h>

SimulationTime* SimulationTime::mpInstance = 0;

/**
 * Return a pointer to the simulation time object.
 * The first time one creates a simulation time object, this method MUST be used.
 * @param durationOfSimulation Total dimensionalized time of the simulation
 * @param totalTimeStepsInSimulation the number of time steps into which the above will be broken
 */
SimulationTime* SimulationTime::Instance(double durationOfSimulation, int totalTimeStepsInSimulation) 
{
	assert(mpInstance == 0);
	mpInstance = new SimulationTime(durationOfSimulation, totalTimeStepsInSimulation);
	return mpInstance;
}

/**
 * Return a pointer to the simulation time object.
 * The second and subsequent times one creates a simulation time object,
 * this method MUST be used.
 */
SimulationTime* SimulationTime::Instance()
{
	assert(mpInstance != 0);
	return mpInstance;
}

SimulationTime::SimulationTime(double durationOfSimulation, int totalTimeStepsInSimulation) {
	mDurationOfSimulation = durationOfSimulation;
	mTotalTimeStepsInSimulation=totalTimeStepsInSimulation;
	mTimeStepsElapsed = 0;
}


void SimulationTime::Destroy() {
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
	return mDurationOfSimulation/mTotalTimeStepsInSimulation;
}

/**
 * Increment the simulation time by one time step
 */
void SimulationTime::IncrementTimeOneStep()
{
	mTimeStepsElapsed++;
}

/**
 * Get the number of time steps that have elapsed.
 * @return number of time steps
 */
int SimulationTime::GetTimeStepsElapsed()
{
	return mTimeStepsElapsed;
}

/**
 * Get the dimensionalised simulation time.
 * Should not have rounding errors.
 * @return simulation time
 */
double SimulationTime::GetDimensionalisedTime()
{
	return ((double)mTimeStepsElapsed / (double)mTotalTimeStepsInSimulation)
	       * mDurationOfSimulation;
}
