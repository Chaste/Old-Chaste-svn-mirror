/**
 * SodiumCurrentLR91.cpp
 * 
 * INa, Fast sodium current.
 */
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
    mMagnitudeOfCurrent = 0.0; // Set to 0.0 initially, it's a variable of IonicCurrent class
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
 */
double SodiumCurrentLR91::GetM()
{   
    return mM;
}

/**
 * Return value of h inactivation gate
 */
double SodiumCurrentLR91::GetH()
{   
    return mH;
}

/**
 * Return value of j inactivation gate
 */
double SodiumCurrentLR91::GetJ()
{   
    return mJ;
}

/**
 * Update gating variables 
 * 
 * @param voltage Current transmembrane voltage
 */
void SodiumCurrentLR91::UpdateAlphaAndBeta(double voltage)
{
    //updating alpha and beta voltage--dependent rate constants
    mAlphaM = 0.32 * (voltage + 47.13)/(1 - exp(-0.1*(voltage + 47.13)));
    mBetaM = 0.08 * exp(-voltage/11.0);
    
    if(voltage >= -40.0)
    {
        mAlphaH = 0.0;
        mBetaH = 1/(0.13*(exp((voltage+10.66)/(-11.1))+1));
        mAlphaJ = 0.0;
        mBetaJ = 0.3*exp(-2.535*0.0000001*voltage)/(1 + exp(-0.1*(voltage+32.0)));
    }
    else
    {
        mAlphaH = 0.135*exp((80+voltage)/(-6.8));
        mBetaH = 3.56*exp(0.079*voltage) + 3.1*0.00001*exp(0.35*voltage);
        mAlphaJ = (-1.2714*0.00001*exp(0.2444*voltage) - 3.474*0.00001*exp(-0.04391*voltage))*(voltage+37.78)/(1+exp(0.311*(voltage+79.23)));
        mBetaJ = 0.1212*exp(-0.0105*voltage)/(1+exp(-0.1378*(voltage+40.14)));
    }
}

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
 */
void SodiumCurrentLR91::UpdateMagnitudeOfCurrent(double voltage, double m, double h, double j)
{   
    SetGatingVariables(m, h, j);
    UpdateAlphaAndBeta(voltage);
    
    // checking how beta is updated
    //std::cout << "\n"<< "mBetaM is equal to " <<mBetaM<<std::endl;
    
    mMagnitudeOfCurrent = gNa * pow(mM, 3) * mH * mJ * (voltage - eNa);
} 

/**
 * Returns m'
 * 
 * @param voltage Current transmembrane potential
 * @param m       Current value for m gating variable
 * @param h       Current value for h gating variable
 * @param j       Current value for j gating variable
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
 */ 
double SodiumCurrentLR91::ComputeJPrime(double voltage, double m, double h, double j)
{       
    SetGatingVariables(m, h, j);
    UpdateAlphaAndBeta(voltage);
    return  (mAlphaJ - (mAlphaJ + mBetaJ)*mJ);
}
