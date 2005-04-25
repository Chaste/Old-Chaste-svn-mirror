#ifndef _BACKGROUNDCURRENTLR91_HPP_
#define _BACKGROUNDCURRENTLR91_HPP_

#include "IonicCurrent.hpp"
#include "ConstantsLR91.hpp"
#include <iostream>

/**
 * Luo-Rudy 91 background current, Ib.
 */
class BackgroundCurrentLR91 : public IonicCurrent
{        
    public:
        // Constructor
        BackgroundCurrentLR91();
        // Destructor
        ~BackgroundCurrentLR91();
        void UpdateMagnitudeOfCurrent(double voltage);
};

#endif //_BACKGROUNDCURRENTLR91_HPP_

