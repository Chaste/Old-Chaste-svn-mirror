/**
 * LR91Model.cpp
 * 
 */
#include "LR91Model.hpp"
#include <cmath>

/**
 * Constructor
 * 
 */
LR91Model::LR91Model(double voltage, double m, double h, double j, double d, 
double f, double x, double caI, AbstractStimulusFunction *pStimulus)
{   
    mV = voltage;
    mM = m;
    mH = h;    
    mJ = j;
    mD = d;
    mF = f;
    mX = x;
    mCaI = caI;
    mpLR91OdeSystem = new LR91OdeFun(pStimulus);
}

/**
 * Destructor
 */
LR91Model::~LR91Model()
{   
    // Do nothing
}


/**
 * Solves the LR91 model using some ODe Solver
 */
std::vector<double> LR91Model::Solve()//tange of time AbstractOdeSolver *pOdeSolver)
{
//    std::vector<double> Y();
//    Y.push_back(mV);
//    Y[1] = mM;
//    Y[2] = mH;
//    Y[3] = mJ;
//    Y[4] = mD;
//    Y[5] = mF;
//    Y[6] = mX;
//    Y[7] = mCaI;
    
   // EvaluateYDerivatives (const double &rTime, std::vector<double> &rY) 
}

