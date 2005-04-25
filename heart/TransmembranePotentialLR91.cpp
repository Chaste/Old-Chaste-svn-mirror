/**
 * TransmembranePotentialLR91.cpp
 * 
 * Voltage
 */
#include "TransmembranePotentialLR91.hpp"
#include <cmath>

/**
 * Constructor
 * 
 */
TransmembranePotentialLR91::TransmembranePotentialLR91()
{   
  //  mV = 0.0;
}

/**
 * Destructor
 */
TransmembranePotentialLR91::~TransmembranePotentialLR91()
{   
    // Do nothing
}

//
//
//void TransmembranePotentialLR91::SetTransmembranePotential(double voltage)
//{
//    mV = voltage;
//}

/**
 * Compute VPrime
 * 
 * @param voltage Current transmembrane voltage
 * @param ITotal  total ionic current
 */
double TransmembranePotentialLR91::ComputeVPrime(double iStim, double iTotal)
{   
          
    return (-(iStim + iTotal)/membraneCapacitance);
} 
