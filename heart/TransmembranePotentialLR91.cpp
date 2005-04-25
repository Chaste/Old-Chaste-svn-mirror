#include "TransmembranePotentialLR91.hpp"
#include <cmath>

/**
 * Constructor
 * 
 */
TransmembranePotentialLR91::TransmembranePotentialLR91()
{   
  //Do nothing
}

/**
 * Destructor
 */
TransmembranePotentialLR91::~TransmembranePotentialLR91()
{   
    // Do nothing
}



/**
 * Compute V'
 * 
 * @param iStim   Stimulus Current
 * @param iTotal  total ionic current
 * 
 * @return double RHS of V'
 */
double TransmembranePotentialLR91::ComputeVPrime(double iStim, double iTotal)
{   
          
    return (-(iStim + iTotal)/membraneCapacitance);
} 
