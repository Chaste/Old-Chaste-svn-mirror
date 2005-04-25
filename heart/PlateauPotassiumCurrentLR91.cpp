#include "PlateauPotassiumCurrentLR91.hpp"
#include <cmath>

/**
 * Constructor
 */
PlateauPotassiumCurrentLR91::PlateauPotassiumCurrentLR91()
{   
    mKp = 0.0;
    mMagnitudeOfCurrent = 0.0;
}

/**
 * Destructor
 */
PlateauPotassiumCurrentLR91::~PlateauPotassiumCurrentLR91()
{   
    // Do nothing
}

/**
 * Update Kp
 * 
 * @param voltage Current transmembrane voltage
 */
void PlateauPotassiumCurrentLR91::UpdateKP(double voltage)
{  
    mKp = 1.0/ ( 1.0 + exp((7.488 - voltage)/5.98) );
}

/**
 * Update magnitude of potassium plateau current, IKp.
 * 
 * @param voltage Current transmembrane voltage
 */
void PlateauPotassiumCurrentLR91::UpdateMagnitudeOfCurrent(double voltage)
{   
    UpdateKP(voltage);        
    mMagnitudeOfCurrent = gKp* mKp * (voltage - eKp);
} 
