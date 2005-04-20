/**
 * IonicConcentrationLR91.cpp
 * 
 * An ionic concentration
 */
#include "IonicConcentrationLR91.hpp"

/**
 * Constructor
 */
IonicConcentrationLR91::IonicConcentrationLR91(void)
{   
    mMagnitudeOfIonicConcentration = 0.0;
}

/**
 * Destructor
 */
IonicConcentrationLR91::~IonicConcentrationLR91(void)
{   
    // Do nothing
}

/**
 * Set magnitude of ionic concentration
 * 
 * @param &rMagnitudeOfCurrent Magnitude of ionic concentration
 */
void IonicConcentrationLR91::SetMagnitudeOfIonicConcentration(const double &rIonicConcentration)
{
    mMagnitudeOfIonicConcentration = rIonicConcentration;
}

/**
 * Get magnitude of ionic concentration 
 */
double IonicConcentrationLR91::GetMagnitudeOfIonicConcentration()
{
    return mMagnitudeOfIonicConcentration;
}
