#include "BackgroundCurrentLR91.hpp"
#include <cmath>

/**
 * Constructor
 */

BackgroundCurrentLR91::BackgroundCurrentLR91()
{   
    mMagnitudeOfCurrent = 0.0; // Set to 0.0 initially, its a variable of IonicCurrent class
}

/**
 * Destructor
 */
BackgroundCurrentLR91::~BackgroundCurrentLR91()
{   
    // Do nothing
}

/**
 * Update magnitude of current, Ib.
 * 
 * @param voltage Current transmembrane voltage
 */
void BackgroundCurrentLR91::UpdateMagnitudeOfCurrent(double voltage)
{   
	// mMagnitudeOfCurrent = gB * (voltage - eB) in original paper, sign of eB changed
	// in ConstantsLR91.hpp, so actually similar
    mMagnitudeOfCurrent = gB * (voltage + eB); 
}
