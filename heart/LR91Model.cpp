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
OdeSolution LR91Model::SolveModel(double startTime, double endTime, double timeStep,
                       std::vector<double> initialConditions, AbstractIvpOdeSolver *pAbstractOdeSolver)
{
     
    OdeSolution  solution = pAbstractOdeSolver->Solve(mpLR91OdeSystem, 
                              startTime, endTime, timeStep, initialConditions);
    return solution;
    
}
                                

