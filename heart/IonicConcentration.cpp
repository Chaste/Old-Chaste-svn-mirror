/**
 * IonicConcentration.cpp
 * 
 * An ionic concentration
 */
#include "IonicConcentration.hpp"

/**
 * Constructor
 */
IonicConcentration::IonicConcentration(void)
{   
    mMagnitudeOfIonicConcentration = 0.0;
}

/**
 * Destructor
 */
IonicConcentration::~IonicConcentration(void)
{   
    // Do nothing
}

/**
 * Set magnitude of ionic concentration
 * 
 * @param &rMagnitudeOfCurrent Magnitude of ionic concentration
 */
void IonicConcentration::SetMagnitudeOfIonicConcentration(const double &rIonicConcentration)
{
    mMagnitudeOfIonicConcentration = rIonicConcentration;
}

/**
 * Get magnitude of ionic concentration 
 */
double IonicConcentration::GetMagnitudeOfIonicConcentration()
{
    return mMagnitudeOfIonicConcentration;
}
