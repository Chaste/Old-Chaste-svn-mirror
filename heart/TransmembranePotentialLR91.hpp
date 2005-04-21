#ifndef _TRANSMEMBRANEPOTENTIALLR91_HPP_
#define _TRANSMEMBRANEPOTENTIALLR91_HPP_


#include "ConstantsLR91.hpp"
#include <iostream>

class TransmembranePotentialLR91 
{
    private:
    // double mV;
               
    public:
        // Constructor
        TransmembranePotentialLR91();
        // Destructor
        ~TransmembranePotentialLR91();
        
        double ComputeVPrime(double iStim, double iTotal);
      //  void SetTransmembranePotential(double voltage);
               
};

#endif //_TRANSMEMBRANEPOTENTIALLR91_HPP_

