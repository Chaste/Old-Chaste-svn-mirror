#ifndef _TRANSMEMBRANEPOTENTIALLR91_HPP_
#define _TRANSMEMBRANEPOTENTIALLR91_HPP_


#include "ConstantsLR91.hpp"
#include <iostream>
/**
 * LR91 Transmembrane Potential V  
 */
class TransmembranePotentialLR91 
{                   
    public:
        // Constructor
        TransmembranePotentialLR91();
        // Destructor
        ~TransmembranePotentialLR91();
        
        double ComputeVPrime(double iStim, double iTotal);
                            
};

#endif //_TRANSMEMBRANEPOTENTIALLR91_HPP_

