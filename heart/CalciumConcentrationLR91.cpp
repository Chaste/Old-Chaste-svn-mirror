#include "CalciumConcentrationLR91.hpp"
#include "SlowInwardCurrentLR91.hpp"
#include <cmath>

/**
 * Constructor
 */
CalciumConcentrationLR91::CalciumConcentrationLR91(void)
{   
    mMagnitudeOfIonicConcentration = 0.0; 
}

/**
 * Destructor
 */
CalciumConcentrationLR91::~CalciumConcentrationLR91()
{   
    // Do nothing
}

/**
 * Compute Ca_i'. In the LR91 model, only intracellular Ca concentration is dynamic.
 * 
 * @param voltage Current transmembrane voltage
 * @param d       Current value of d gating variable
 * @param f       Current value of f gating variable
 * @param caI     Intracellular calcium concentration
 * @param iSi     Current value of slow inward current
 * 
 * @return double The updated intracellular Ca concentration
 */
double CalciumConcentrationLR91::ComputeCalciumPrime(double voltage, double d, double f, double caI, double iSi)
{   
	SetMagnitudeOfIonicConcentration(caI);
    return (-1e-4 * iSi + 0.07 * (1e-4 -  mMagnitudeOfIonicConcentration));   
}
