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
LR91Model::LR91Model(double voltage, double m, double h, double j, double d, double f, double x, double caI)
{   
    mV = voltage;
    mM = m;
    mH = h;    
    mJ = j;
    mD = d;
    mF = f;
    mX = x;
    mCaI = caI;
    mpLR91OdeSystem = new LR91OdeFun();
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
std::vector<double> LR91Model::Solve()//OdeSolver myOdeSolver)
{
}

