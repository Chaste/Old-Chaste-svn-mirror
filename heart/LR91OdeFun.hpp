#ifndef _LR91ODEFUN_HPP_
#define _LR91ODEFUN_HPP_

#include "ConstantsLR91.hpp"

//#include "AbstractOdeSystem.hpp"

#include "TransmembranePotentialLR91.hpp"
#include "SodiumCurrent.hpp"
//#include "SodiumCurrentLR91.hpp"
#include "BackgroundCurrentLR91.hpp"
#include "CalciumConcentrationLR91.hpp"
#include "PlateauPotassiumCurrentLR91.hpp"
#include "SlowInwardCurrentLR91.hpp"
//#include "PotassiumTimeInpendententCurrentLR91.hpp"
//#include "PotassiumTimeDependentCurrentLR91.hpp"
 
#include <iostream>

//#include "petscvec.h"

class LR91OdeFun //: public AbstractOdeSystem
{
    private:
        TransmembranePotentialLR91 *mpV;
        SodiumCurrent  *mpINa;
        //SodiumCurrentLR91  *mpINa;
        BackgroundCurrentLR91 *mpIB;
        CalciumConcentrationLR91 *mpCaI;
        PlateauPotassiumCurrentLR91 *mpKP;
        SlowInwardCurrentLR91 *mpISi;
        //PotassiumTimeInpendententCurrentLR91 *mpIK1;
        //PotassiumTimeDependentCurrentLR91 *mpIK;
       
          
    public:
        // Constructor
        LR91OdeFun();
        // Destructor
        ~LR91OdeFun();
        
        // This mehtod will compute the RHS of the Luo--Rudy model
        // void EvaluateYDiffs(double rTime, Vec rY, Vec rYNew);   
        //void ComputingRHS(tOfStimulus, std::vector<double> initCond); 
};
 
#endif //_LR91ODEFUN_HPP_
