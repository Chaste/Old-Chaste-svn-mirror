/**
 * CalciumConcentrationLR91.cpp
 * 
 * [Ca_i], Intracellular Calcium concentration
 */
#include "CalciumConcentrationLR91.hpp"
#include "SlowInwardCurrentLR91.hpp"
#include <cmath>

/**
 * Constructor
 */
CalciumConcentrationLR91::CalciumConcentrationLR91(void)
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
 * Compute Ca_i'.
 * 
 * @param voltage Current transmembrane voltage
 * @param d       Current value of d gating variable
 * @param f       Current value of f gating variable
 * @param caI     Intracellular calcium concentration
 * @param iSi     Current value of slow inward current
 */
double CalciumConcentrationLR91::ComputeCalciumPrime(double voltage, double d, double f, double caI, double iSi)
{   
//    SlowInwardCurrentLR91 *pISi;
//    pISi = new SlowInwardCurrentLR91();
//    pISi->UpdateMagnitudeOfCurrent(voltage,d,f,caI);
//    double iSi = pISi->GetMagnitudeOfCurrent();
        
    SetMagnitudeOfIonicConcentration(caI);
    return (-0.0001 * iSi + 0.007 * (0.0001 -  mMagnitudeOfIonicConcentration));   
}
