#ifndef _LR91ODEFUN_HPP_
#define _LR91ODEFUN_HPP_

#include "ConstantsLR91.hpp"
//#include "AbstractOdeSystem.hpp"
#include "SodiumCurrent.hpp"
#include "BackgroundCurrentLR91.hpp"
#include <iostream>

#include "petscvec.h"

class LR91OdeFun //: public AbstractOdeSystem
{
    public:
        // Constructor
        LR91OdeFun();
        // Destructor
        ~LR91OdeFun();
        
        // This mehtod will compute the RHS of the Luo--Rudy model
        void EvaluateYDiffs(double rTime, Vec rY, Vec rYNew);   
};
 
#endif //_LR91ODEFUN_HPP_
