/**
 * SlowInwardCurrentLR91.cpp
 * 
 * ISi.
 */
#include "SlowInwardCurrentLR91.hpp"
#include <cmath>

/**
 * Constructor
 */
SlowInwardCurrentLR91::SlowInwardCurrentLR91()
{   
    mD = 0.0;
    mF = 0.0;
    mESi = 0.0;    
    mAlphaD = 0.0;
    mAlphaF = 0.0;
    mBetaD = 0.0;
    mBetaF = 0.0;
    mMagnitudeOfCurrent = 0.0; // Set to 0.0 initially, its a variable of IonicCurrent class
}

/**
 * Destructor
 */
SlowInwardCurrentLR91::~SlowInwardCurrentLR91()
{   
    // Do nothing
}

/**
 * Return value of d gate
 */
double SlowInwardCurrentLR91::GetD()
{   
    return mD;
}

/**
 * Return value of f gate
 */
double SlowInwardCurrentLR91::GetF()
{   
    return mF;
}



/**
 * Update voltage-dependent rate constants, alpha and beta
 * 
 * @param voltage Current transmembrane voltage
 */
void SlowInwardCurrentLR91::UpdateAlphaAndBeta(double voltage)
{
    double mAlphaD = 0.095 * exp(-0.01*(voltage - 5.0)) / (1 + exp(-0.072*(voltage - 5.0)));
    double mBetaD =   0.07 * exp(-0.017*(voltage + 44.0)) / (1 + exp(-0.05*(voltage + 44.0)));
    double mAlphaF = 0.012 * exp(-0.008*(voltage + 28.0)) / (1 + exp(-0.15*(voltage + 28.0)));
    double mBetaF = 0.0065 * exp(0.02 * (voltage + 30.0)) /  (1 + exp(-0.2*(voltage + 30.0)));
}


/**
 * Set gating variables to given values
 * 
 * @param d New value for d gating variable
 * @paran f New value for f gating variable
 */
void SlowInwardCurrentLR91::SetGatingVariables(double d, double f)
{  
    mD = d;
    mF = f;
}

/**
 * Update ESi, slow inward current reversal potential
 * 
 * @param caI [Ca_i] concentration
 */
void SlowInwardCurrentLR91::UpdateESi(double caI)
{
    mESi = 7.7 - 13.0287 * log(caI);    
}

/**
 * Update magnitude of current, ISi.
 * 
 * @param voltage Current transmembrane voltage
 */
void SlowInwardCurrentLR91::UpdateMagnitudeOfCurrent(double voltage, double d, double f, double caI)
{   
    SetGatingVariables(d, f);
    UpdateAlphaAndBeta(voltage);
    UpdateESi(caI);    
    mMagnitudeOfCurrent = gSi * mD * mF * (voltage - mESi);
} 

/**
 * Returns y'
 * 
 * @param voltage Current transmembrane potential
 * @param d       Current value for d gating variable
 * @param f       Current value for f gating variable
 */
double SlowInwardCurrentLR91::ComputeDPrime(double voltage, double d, double f)
{    
    SetGatingVariables(d, f);
    UpdateAlphaAndBeta(voltage);
    return (mAlphaD - (mAlphaD + mBetaD)*mD); 
}

/**
 * Returns f'
 * 
 * @param voltage Current transmembrane potential
 * @param d       Current value for d gating variable
 * @param f       Current value for f gating variable
 */
double SlowInwardCurrentLR91:: ComputeFPrime(double voltage, double d, double f)
{    
    SetGatingVariables(d, f);
    UpdateAlphaAndBeta(voltage);
    return  (mAlphaF - (mAlphaF + mBetaF)*mF);
}
