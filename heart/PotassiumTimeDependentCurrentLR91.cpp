/**
 * PotassiumTimeDependentCurrentLR91.cpp
 * 
 * Ik, Time dependent potassium current.
 */
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
    mMagnitudeOfCurrent = 0.0; // Set to 0.0 initially, its a variable of IonicCurrent class
}

/**
 * Destructor
 */
PotassiumTimeDependentCurrentLR91::~PotassiumTimeDependentCurrentLR91()
{   
    // Do nothing
}

/**
 * Return value of x activation gate
 */
double PotassiumTimeDependentCurrentLR91::GetX()
{   
    return mX;
}

double PotassiumTimeDependentCurrentLR91::GetXi(double voltage)
{
	UpdateXi(voltage) ;
	return mXi;	
}


/**
 * Update gating variables 
 * 
 * @param voltage Current transmembrane voltage
 */
void PotassiumTimeDependentCurrentLR91::UpdateAlphaAndBeta(double voltage)
{
    //updating alpha and beta voltage--dependent rate constants
    double mAlphaX= 0.0005*exp(0.083*(voltage + 50.0))/ (1 + exp(0.057*(voltage + 50.0))) ;
    double mBetaX =  0.0013*exp(-0.06*(voltage +20.0))/(1 +exp(-0.04*(voltage + 20.0)) ) ;
}

/**
 * Update Xi
 * 
 * @param voltage Current transmembrane potential
 */
void PotassiumTimeDependentCurrentLR91::UpdateXi(double voltage)
{
    //updating Xi voltage--dependent rate constants
    
    if (voltage <= -100.0)
    {
    	double mXi= 1.0  ;
    }
    else
    {
     double mXi = 2.837*(exp(0.04*(voltage+77.0)) -1) / ( (voltage +77.0)*exp(0.04*(voltage+35.0)) ) ;	
    } 
}

/**
 * Set gating variables to specified values
 * 
 * @param x New value for X gating variable
 */
void PotassiumTimeDependentCurrentLR91::SetGatingVariables(double x)
{  
    mX = x; 
}

/**
 * Update magnitude of current, IK.
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
 * @param voltage Current transmembrane potential
 * @param x       Current value for X gating variable
 */
double PotassiumTimeDependentCurrentLR91::ComputeXPrime(double voltage, double x)
{    
    SetGatingVariables(x);
    UpdateAlphaAndBeta(voltage);
    UpdateXi(voltage);
    return  (mAlphaX - (mAlphaX + mBetaX)*mX); 
}
