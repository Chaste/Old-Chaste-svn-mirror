/**
 * PotassiumTimeIndependentCurrentLR91.cpp
 * 
 * IK1, Time independent potassium current.
 */
#include "PotassiumTimeIndependentCurrentLR91.hpp"
#include <cmath>

/**
 * Constructor
 */
PotassiumTimeIndependentCurrentLR91::PotassiumTimeIndependentCurrentLR91()
{   
    mK1 = 0.0; 
    mAlphaK1 = 0.0;
    mBetaK1 = 0.0;
    mMagnitudeOfCurrent = 0.0; // Set to 0.0 initially, its a variable of IonicCurrent class
}

/**
 * Destructor
 */
PotassiumTimeIndependentCurrentLR91::~PotassiumTimeIndependentCurrentLR91()
{   
    // Do nothing
}

/**
 * Update alpha and beta variables 
 * 
 * @param voltage Current transmembrane voltage
 */
void PotassiumTimeIndependentCurrentLR91::UpdateAlphaAndBeta(double voltage)
{
    //updating alpha and beta voltage-dependent rate constants
    mAlphaK1 = 1.02 / ( 1.0 + exp(0.2385*(voltage - eK1 - 59.215 )) ) ;
    mBetaK1 =  (0.49124*exp(0.08032*(voltage-eK1 + 5.476)) + exp(0.06175*(voltage - eK1 - 594.31)))/ (1+exp(-0.5143*(voltage-eK1+4.753))) ;   
}

/**
 * Update magnitude of current, IK1.
 * 
 * @param voltage Current transmembrane voltage
 */
void PotassiumTimeIndependentCurrentLR91::UpdateMagnitudeOfCurrent(double voltage)
{   
    // UpdateAlphaAndBeta(voltage);
    UpdateK1(voltage);
    mMagnitudeOfCurrent = gK1*mK1*(voltage - eK1);
}

/**
 * Update K1 inactivation gate
 * 
 * @param voltage Current transmembrane voltage
 */
void PotassiumTimeIndependentCurrentLR91::UpdateK1(double voltage)
{
	UpdateAlphaAndBeta(voltage);
	mK1 = (mAlphaK1)/(mAlphaK1 + mBetaK1);
}

double PotassiumTimeIndependentCurrentLR91::GetK1(double voltage)
{
	UpdateAlphaAndBeta(voltage) ;
	return mK1 ;
}
