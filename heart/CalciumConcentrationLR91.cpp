/**
 * CalciumConcentrationLR91.cpp
 * 
 * Cai.
 */
#include "CalciumConcentrationLR91.hpp"
#include "SlowInwardCurrentLR91.hpp"
#include <cmath>

/**
 * Constructor
 * 
 * @param CaI    Initial value for intracellular calcium concentration
 */
CalciumConcentrationLR91::CalciumConcentrationLR91()
{   
    mMagnitudeOfIonicConcentration = 0.0; // Set to 0.0 initially, its a variable of IonicCurrent class
}

/**
 * Destructor
 */
CalciumConcentrationLR91::~CalciumConcentrationLR91()
{   
    // Do nothing
}


/**
 * Get Kp.
 * 
 */

double CalciumConcentrationLR91::GetCaI()
{   
    return mMagnitudeOfIonicConcentration;
}


/**
 * Update magnitude of intracellular Ca concentration, Cai.
 * 
 * @param voltage Current transmembrane voltage
 * @param caI     Intracellular calcium concentration
 */
void CalciumConcentrationLR91::UpdateMagnitudeOfIonicConcentration(double caI)
{    
    mMagnitudeOfIonicConcentration = caI;
} 

/**
 * Compute CaiPrime.
 * 
 * @param voltage Current transmembrane voltage
 * @param caI     Intracellular calcium concentration
 */
double CalciumConcentrationLR91::ComputeCalciumPrime(double voltage, double d, double f, double caI, double iSi)
{   
//    SlowInwardCurrentLR91 *pISi;
//    pISi = new SlowInwardCurrentLR91();
//    pISi->UpdateMagnitudeOfCurrent(voltage,d,f,caI);
//    double iSi = pISi->GetMagnitudeOfCurrent();
        
    UpdateMagnitudeOfIonicConcentration(caI);
    return (-0.0001 * iSi + 0.007 * (0.0001 -  mMagnitudeOfIonicConcentration));   
}
