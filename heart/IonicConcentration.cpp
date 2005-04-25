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
 * 
 * @return double Magnitude of Ionic Concentration
 */
double IonicConcentration::GetMagnitudeOfIonicConcentration(void)
{
    return mMagnitudeOfIonicConcentration;
}
