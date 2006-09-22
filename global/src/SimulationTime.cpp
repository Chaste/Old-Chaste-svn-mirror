/**
 * Simulation time object stores the simulation time, uses the singleton pattern
 */

#include "SimulationTime.hpp"

SimulationTime* SimulationTime::mInstance = 0;

/**
 * Return a pointer to the simulation time object
 */
SimulationTime* SimulationTime::Instance() 
{
	if (mInstance == 0)
	{ 
		mInstance = new SimulationTime();
	}
	return mInstance;
}

SimulationTime::SimulationTime() {
	mTime=0;
}

/**
 * Get the simlation time
 * @return simluation time
 */
double SimulationTime::GetTime()
{
	return mTime;
}
/**
 * Increment the simulation time
 * @param dt the amount added to the time
 */
void SimulationTime::IncrementTime(double dt)
{
	mTime+=dt;
}

