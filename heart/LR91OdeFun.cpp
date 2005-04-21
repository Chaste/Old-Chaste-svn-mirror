/**
 * LR91OdeFun.cpp
 * 
 */
#include "LR91OdeFun.hpp"
#include <cmath>

/**
 * Constructor
 * 
 */
LR91OdeFun::LR91OdeFun()
{
    mpV = new TransmembranePotentialLR91();
    mpINa = new SodiumCurrent();
    //mpINa = new SodiumCurrentLR91();
    mpIB = new BackgroundCurrentLR91();
    mpCaI = new CalciumConcentrationLR91();
    mpIKp = new PlateauPotassiumCurrentLR91();
    mpISi = new SlowInwardCurrentLR91();
    //mIK1 = new PotassiumTimeInpendententCurrentLR91();;
    //mIK =  new SPotassiumTimeDependentCurrentLR91();  
    
}

/**
 * Destructor
 */
LR91OdeFun::~LR91OdeFun()
{   
    // Do nothing
}

/**
 * Method that evaluates the RHS of the Luo--Rudy model 
 *  ??! are rYNew updated within, i.e. what it points to is modified?! 
 * 
// */
//void LR91OdeFun::EvaluateYDiffs(double rTime, Vec rY, Vec rYNew)
//{
////    // need to extract V, m,h,...from Petski vector
////   Vec mX = rY;
////   
////   PetscScalar *xElements;
////
////   VecGetArray(mX,xElements);
////   
////   xElements[0] = rY[0];
////   double V = xElements[0];
////   
////   double m =xElements[1];
////   
////   
////    
//    
//    // update all currents, gate variables, cocentrations 
//    
//    // compute RHS and put it into Petski vector form
//       
////}
//
//void LR91OdeFun::ComputingRHS(double tOfStimulus, std::vector<double> initCond)
std::vector<double> LR91OdeFun::EvaluateYDerivatives (const double &rTime, std::vector<double> &rY) 
{
//    double v = initCond[0];
//    double m = initCond[1];
//    double h = initCond[2];
//    double j = initCond[3];
//    double d = initCond[4];
//    double f = initCond[5];
//    double x = initCond[6];
//    double caI = initCond[7];   
    double v = rY[0];
    double m = rY[1];
    double h = rY[2];
    double j = rY[3];
    double d = rY[4];
    double f = rY[5];
    double x = rY[6];
    double caI = rY[7];
    

    
    // Compute all the currents
    //sodium current
    mpINa->UpdateMagnitudeOfCurrent(v, m, h, j);
    double iNa = mpINa->GetMagnitudeOfCurrent();
    
    //slow inward current
    mpISi->UpdateMagnitudeOfCurrent(v,d,f,caI);
    double iSi = mpISi->GetMagnitudeOfCurrent();
    
    //potassium time dependent current
    //mpIK->UpdateMagnitudeOfCurrent(v);
    //double iK = mpIK->GetMagnitudeOfCurrent();
    
    //potassium time independent current
    //mpIK1->UpdateMagnitudeOfCurrent(v);
    //double iK1 = mpIK1->GetMagnitudeOfCurrent();
    
    //potassium plateau
    mpIKp->UpdateMagnitudeOfCurrent(v);
    double iKp = mpIKp->GetMagnitudeOfCurrent();
    
    //bacground current
    mpIB->UpdateMagnitudeOfCurrent(v);
    double iB = mpIB->GetMagnitudeOfCurrent();
    
    //Update Ca 
    mpCaI->SetMagnitudeOfIonicConcentration(caI);
    double iCa = mpCaI->GetMagnitudeOfIonicConcentration();
        
    
    // Compute All the RHSs
    
    //RHS of gating variables
    double mPrime = mpINa->ComputeMPrime(v, m, h, j);
    double hPrime = mpINa-> ComputeHPrime(v, m, h, j);
    double jPrime = mpINa->ComputeJPrime(v, m, h, j);
    double dPrime = mpISi->ComputeDPrime(v, d, f);
    double fPrime = mpISi->ComputeFPrime(v, d, f);
    double xPrime = 0.0;//pIK->ComputeFPrime(v);
    double caIPrime = mpCaI->ComputeCalciumPrime(v, d,f, caI, iSi);
        
    // Total Current
    double iTotal = iSi + iNa + iKp + iB; //+ iK + iK1;
    
    // Introducing Stimulus -- need to modify
    double iStim = 0.0;
    if (rTime==1)
    {
        iStim = 20.0;
    }
    else 
    {
        iStim = 0.0;
    }
    
    // Calculating VPrime
    double VPrime = mpV->ComputeVPrime(iStim, iTotal);
    
    std::vector<double> returnRHS;
    returnRHS[0] = VPrime;
    returnRHS[1] = mPrime;
    returnRHS[2] = hPrime;
    returnRHS[3] = jPrime;
    returnRHS[4] = dPrime;
    returnRHS[5] = fPrime;
    returnRHS[6] = xPrime;
    returnRHS[7] = caIPrime;
    
    return returnRHS;
}

