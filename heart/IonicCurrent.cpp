#include "IonicCurrent.hpp"

/**
 * Constructor
 */
IonicCurrent::IonicCurrent(void)
{	
	mMagnitudeOfCurrent = 0.0;
}

/**
 * Overloaded Constructor
 * 
 * @param &rMagnitudeOfCurrent Initial magnitude of current
 */
IonicCurrent::IonicCurrent(const double &rMagnitudeOfCurrent)
{	
	mMagnitudeOfCurrent = rMagnitudeOfCurrent;
}

/**
 * Destructor
 */
IonicCurrent::~IonicCurrent(void)
{	
	// Do nothing
}

/**
 * Set magnitude of current
 * 
 * @param &rMagnitudeOfCurrent Magnitude of current
 */
void IonicCurrent::SetMagnitudeOfCurrent(const double &rMagnitudeOfCurrent)
{
	mMagnitudeOfCurrent = rMagnitudeOfCurrent;
}

/**
 * Get magnitude of current
 * 
 * @return double Magnitude of current
 */
double IonicCurrent::GetMagnitudeOfCurrent(void)
{
	return mMagnitudeOfCurrent;
}
