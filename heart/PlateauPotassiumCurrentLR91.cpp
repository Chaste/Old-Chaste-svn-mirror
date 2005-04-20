/**
 * PlateauPotassiumCurrentLR91.cpp
 * 
 * IKp.
 */
#include "PlateauPotassiumCurrentLR91.hpp"
#include <cmath>

/**
 * Constructor
 * 
 * @param Kp Initial value for Kp gate
 */
PlateauPotassiumCurrentLR91::PlateauPotassiumCurrentLR91()
{   
    mKp = 0.0;
    mMagnitudeOfCurrent = 0.0; // Set to 0.0 initially, its a variable of IonicCurrent class
}

/**
 * Destructor
 */
PlateauPotassiumCurrentLR91::~PlateauPotassiumCurrentLR91()
{   
    // Do nothing
}


/**
 * Get Kp.
 * 
 */

double PlateauPotassiumCurrentLR91::GetKp()
{   
    return mKp;
}


/**
 * Update Kp
 * 
 * @param voltage Current transmembrane voltage
 */

void PlateauPotassiumCurrentLR91::UpdateGatingVariables(double voltage)
{  
    mKp = 1/ ( 1 + exp(7.488 - voltage)/5.98);
}


/**
 * Update magnitude of current, IKp.
 * 
 * @param voltage Current transmembrane voltage
 */
void PlateauPotassiumCurrentLR91::UpdateMagnitudeOfCurrent(double voltage)
{   
    UpdateGatingVariables(voltage);        
    mMagnitudeOfCurrent = gKp* mKp * (voltage - eKp);
} 


