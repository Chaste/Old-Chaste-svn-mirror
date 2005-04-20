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
    // nothing to initialize ?!? since no private variables !?!...why then have mYInit... (check with AbstractOdeSystem.hpp) 
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
 */
void EvaluateYDiffs(double rTime, Vec rY, Vec rYNew)
{
    // need to extract V, m,h,...from Petski vector
    
    // update all currents, gate variables, cocentrations 
    
    // compute RHS and put it into Petski vector form
       
}

