#include "PotassiumTimeDependentCurrentLR91.hpp"
#include <cmath>

/**
 * Constructor
 */
PotassiumTimeDependentCurrentLR91::PotassiumTimeDependentCurrentLR91()
{   
    mX = 0.0;
    mAlphaX = 0.0;
    mBetaX = 0.0;
    mXi = 0.0 ;
    mMagnitudeOfCurrent = 0.0;
}

/**
 * Destructor
 */
PotassiumTimeDependentCurrentLR91::~PotassiumTimeDependentCurrentLR91()
{   
    // Do nothing
}

/**
 * Get value of X activation gate
 * 
 * @return double Value of X activation gate
 */
double PotassiumTimeDependentCurrentLR91::GetX()
{   
    return mX;
}

/**
 * Get value of Xi, voltage dependent term of Ik
 * 
 * @return double Value of Xi
 */
double PotassiumTimeDependentCurrentLR91::GetXi(double voltage)
{
	UpdateXi(voltage) ;
	return mXi;	
}

/**
 * Update alpha and beta voltage-dependent rate variables.
 * 
 * @param voltage Current transmembrane voltage
 */
void PotassiumTimeDependentCurrentLR91::UpdateAlphaAndBeta(double voltage)
{
    mAlphaX= 0.0005*exp(0.083*(voltage + 50.0))/ (1 + exp(0.057*(voltage + 50.0))) ;
    mBetaX =  0.0013*exp(-0.06*(voltage + 20.0))/(1 +exp(-0.04*(voltage + 20.0)) ) ;
}

/**
 * Update Xi, voltage dependent term of Ik
 * 
 * @param voltage Current transmembrane voltage
 */
void PotassiumTimeDependentCurrentLR91::UpdateXi(double voltage)
{   
    if (voltage <= -100.0)
    {
    	mXi= 1.0;
    }
    else
    {
      mXi = 2.837*(exp(0.04*(voltage+77.0)) - 1) / ( (voltage +77.0)*exp(0.04*(voltage+35.0)) );	
    } 
}

/**
 * Set gating variable, X
 * 
 * @param x New value for X gating variable
 */
void PotassiumTimeDependentCurrentLR91::SetGatingVariables(double x)
{  
    mX = x; 
}

/**
 * Update magnitude of current, Ik.
 * 
 * @param voltage Current transmembrane voltage
 * @param x       Current value of X gating variable
 */
void PotassiumTimeDependentCurrentLR91::UpdateMagnitudeOfCurrent(double voltage, double x)
{   
    SetGatingVariables(x);
    UpdateAlphaAndBeta(voltage);
    UpdateXi(voltage);
    mMagnitudeOfCurrent = gK*mX*mXi*(voltage - eK);
} 

/**
 * Returns X'
 * 
 * @param voltage Current transmembrane voltage
 * @param x       Current value for X gating variable
 * 
 * @return double RHS of X'
 */
double PotassiumTimeDependentCurrentLR91::ComputeXPrime(double voltage, double x)
{    
    SetGatingVariables(x);
    UpdateAlphaAndBeta(voltage);
    return (mAlphaX - (mAlphaX + mBetaX)*mX); 
}
