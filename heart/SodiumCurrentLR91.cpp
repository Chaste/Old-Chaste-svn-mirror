#include "SodiumCurrentLR91.hpp"
#include <cmath>

/**
 * Constructor
 * 
 * @param m Initial value for m activation gate
 * @param h Initial value for h inactivation gate
 * @param j Initial value for j inactivation gate
 */
SodiumCurrentLR91::SodiumCurrentLR91()
{   
    mM = 0.0;
    mH = 0.0;
    mJ = 0.0;
    
    mAlphaM = 0.0;
    mAlphaH = 0.0;
    mAlphaJ = 0.0;
    mBetaM = 0.0;
    mBetaH = 0.0;
    mBetaJ = 0.0;
    mMagnitudeOfCurrent = 0.0; 
}

/**
 * Destructor
 */
SodiumCurrentLR91::~SodiumCurrentLR91()
{   
    // Do nothing
}

/**
 * Return value of m activation gate
 * 
 * @return double value of m gating variable
 */
double SodiumCurrentLR91::GetM()
{   
    return mM;
}

/**
 * Return value of h inactivation gate
 * 
 * @return double value of h gating variable
 */
double SodiumCurrentLR91::GetH()
{   
    return mH;
}

/**
 * Return value of j inactivation gate
 * 
 * @return  double value of j gating variable
 */
double SodiumCurrentLR91::GetJ()
{   
    return mJ;
}

/**
 * Update alpha and beta voltage-dependent rate variables 
 * 
 * @param voltage Current transmembrane voltage
 */
void SodiumCurrentLR91::UpdateAlphaAndBeta(double voltage)
{
    mAlphaM = 0.32 * (voltage + 47.13)/(1 - exp(-0.1*(voltage + 47.13)));
    mBetaM = 0.08 * exp(-voltage/11.0);
    
    if(voltage >= -40.0)
    {
        mAlphaH = 0.0;
        mBetaH = 1 / (0.13*(1.0+exp((voltage+10.66)/(-11.1))));
        mAlphaJ = 0.0;
        mBetaJ = 0.3*exp(-2.535e-7*voltage)/(1.0 + exp(-0.1*(voltage+32.0)));
    }
    else
    {
        mAlphaH = 0.135*exp((80+voltage)/(-6.8));
        mBetaH = 3.56*exp(0.079*voltage) + 3.1e5*exp(0.35*voltage);
        mAlphaJ = (-1.2714e5*exp(0.2444*voltage) - 3.474e-5*exp(-0.04391*voltage))*(voltage+37.78)/(1+exp(0.311*(voltage+79.23)));
        mBetaJ = 0.1212*exp(-0.01052*voltage)/(1+exp(-0.1378*(voltage+40.14)));
    }
}

/**
 * Set Gating variables
 * 
 * @param m  gating variable m 
 * @param h  gating variable h
 * @param j  gating variable j 
 */

void SodiumCurrentLR91::SetGatingVariables(double m, double h, double j)
{  
    mM = m;
    mH = h;
    mJ = j;
}

/**
 * Update magnitude of current, INa.
 * 
 * @param voltage Current transmembrane voltage
 * @param m  gating variable m 
 * @param h  gating variable h
 * @param j  gating variable j
 */
void SodiumCurrentLR91::UpdateMagnitudeOfCurrent(double voltage, double m, double h, double j)
{   
    SetGatingVariables(m, h, j);         
    mMagnitudeOfCurrent = gNa * pow(mM, 3) * mH * mJ * (voltage - eNa);
} 

/**
 * Returns m'
 * 
 * @param voltage Current transmembrane potential
 * @param m       Current value for m gating variable
 * @param h       Current value for h gating variable
 * @param j       Current value for j gating variable
 * 
 * @return double  rhs of m' ODE 
 */ 
double SodiumCurrentLR91::ComputeMPrime(double voltage, double m, double h, double j)
{    
    SetGatingVariables(m, h, j);
    UpdateAlphaAndBeta(voltage);
    return  (mAlphaM - (mAlphaM + mBetaM)*mM); 
}

/**
 * Returns h'
 * 
 * @param voltage Current transmembrane potential
 * @param m       Current value for m gating variable
 * @param h       Current value for h gating variable
 * @param j       Current value for j gating variable
 * 
 * @return double rhs of h' ODE
 */ 
double SodiumCurrentLR91::ComputeHPrime(double voltage, double m, double h, double j)
{    
    SetGatingVariables(m, h, j);
    UpdateAlphaAndBeta(voltage);
    return  (mAlphaH - (mAlphaH + mBetaH)*mH);
}

/**
 * Returns j'
 * 
 * @param voltage Current transmembrane potential
 * @param m       Current value for m gating variable
 * @param h       Current value for h gating variable
 * @param j       Current value for j gating variable
 * 
 * @return double rhs of j' ODE
 */ 
double SodiumCurrentLR91::ComputeJPrime(double voltage, double m, double h, double j)
{       
    SetGatingVariables(m, h, j);
    UpdateAlphaAndBeta(voltage);
    return  (mAlphaJ - (mAlphaJ + mBetaJ)*mJ);
}
