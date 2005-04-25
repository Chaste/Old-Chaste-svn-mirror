/**
 * LR91OdeFun.cpp
 * 
 */
#include "LR91OdeFun.hpp"
#include <cmath>
#include <cassert>
/**
 * Constructor
 * 
 */
LR91OdeFun::LR91OdeFun(AbstractStimulusFunction *stimulus): AbstractOdeSystem(8)
{
    mpV = new TransmembranePotentialLR91();
    mpINa = new SodiumCurrentLR91();
    mpIB = new BackgroundCurrentLR91();
    mpCaI = new CalciumConcentrationLR91();
    mpIKp = new PlateauPotassiumCurrentLR91();
    mpISi = new SlowInwardCurrentLR91();
    mpIK1 = new PotassiumTimeIndependentCurrentLR91();
    mpIK =  new PotassiumTimeDependentCurrentLR91();  
    mpStimulus = stimulus;
}

/**
 * Destructor
 */
LR91OdeFun::~LR91OdeFun()
{   
    // Do nothing
}

std::vector<double> LR91OdeFun::EvaluateYDerivatives (double time, const std::vector<double> &rY) 
{  
     /*
     * Throw an exception if the initial vector is larger than size 8
     * 
     */
    
    assert(rY.size() == 8);
  
    double v = rY[0];
    double m = rY[1];
    double h = rY[2];
    double j = rY[3];
    double d = rY[4];
    double f = rY[5];
    double x = rY[6];
    double caI = rY[7];
    
    /*
     * Assert that gating variables are within [0,1] range
     */
    assert(m <= 1 && m >= 0 && h <= 1 &&  h >= 0 && j <= 1 &&  j >= 0 && d <= 1 &&  d >= 0 && f <= 1 &&  f >= 0);
 
    // Compute all the currents
    //sodium current
    mpINa->UpdateMagnitudeOfCurrent(v, m, h, j);
    double iNa = mpINa->GetMagnitudeOfCurrent();
   
    //slow inward current
    mpISi->UpdateMagnitudeOfCurrent(v,d,f,caI);
    double iSi = mpISi->GetMagnitudeOfCurrent();
     
    //potassium time dependent current
    mpIK->UpdateMagnitudeOfCurrent(v, x);
    double iK = mpIK->GetMagnitudeOfCurrent();
    
    //potassium time independent current
    mpIK1->UpdateMagnitudeOfCurrent(v);
    double iK1 = mpIK1->GetMagnitudeOfCurrent();
    
    //potassium plateau
    mpIKp->UpdateMagnitudeOfCurrent(v);
    double iKp = mpIKp->GetMagnitudeOfCurrent();
    
    //bacground current
    mpIB->UpdateMagnitudeOfCurrent(v);
    double iB = mpIB->GetMagnitudeOfCurrent();
    
    // Compute All the RHSs
    
    //RHS of gating variables
    double mPrime = mpINa->ComputeMPrime(v, m, h, j);
    double hPrime = mpINa-> ComputeHPrime(v, m, h, j);
    double jPrime = mpINa->ComputeJPrime(v, m, h, j);
    double dPrime = mpISi->ComputeDPrime(v, d, f);
    double fPrime = mpISi->ComputeFPrime(v, d, f);
    double xPrime = mpIK->ComputeXPrime(v, x); 
    double caIPrime = mpCaI->ComputeCalciumPrime(v, d,f, caI, iSi);
        
    // Total Current
    double iTotal = iSi + iNa + iKp + iB + iK + iK1;
    
    double iStim = mpStimulus->GetStimulus(time);
    
    // Calculating VPrime
    double VPrime = mpV->ComputeVPrime(iStim, iTotal); 
    
  
    
    std::vector<double> returnRHS;
    returnRHS.push_back(VPrime);
    returnRHS.push_back(mPrime);
    returnRHS.push_back(hPrime);
    returnRHS.push_back(jPrime);
    returnRHS.push_back(dPrime);
    returnRHS.push_back(fPrime);
    returnRHS.push_back(xPrime);
    returnRHS.push_back(caIPrime);
    

    
    return returnRHS;
}

