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
    mpKP = new PlateauPotassiumCurrentLR91();
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
void LR91OdeFun::ComputingRHS(double tOfStimulus, std::vector<double> initCond)
{
    double v = initCond[0];
    double m = initCond[1];
    double h = initCond[2];
    double j = initCond[3];
    double d = initCond[4];
    double f = initCond[5];
    double x = initCond[6];
    double caI = initCond[7];
    
    mpISi->UpdateMagnitudeOfCurrent(v,d,f,caI);
    double iSi = mpISi->GetMagnitudeOfCurrent();
    
    double iTotal = iSi; //+ iNa + iK + iK1 + iKP + iB;
    double iStim = 0.0;
    if (tOfStimulus==1)
    {
        iStim = 20.0;
    }
    else 
    {
        iStim = 0.0;
    }
    
        double VPrime = mpV->ComputeVPrime(iStim, iTotal);
        
}

