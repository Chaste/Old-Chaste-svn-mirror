#ifndef _LR91ODEFUN_HPP_
#define _LR91ODEFUN_HPP_

#include "ConstantsLR91.hpp"
//#include "AbstractOdeSystem.hpp"
#include "SodiumCurrentLR91.hpp"
#include "BackgroundCurrentLR91.hpp"
#include <iostream>

#include "petscvec.h"

class LR91OdeFun //: public AbstractOdeSystem
{
<<<<<<< .mine
=======
    private:
        TransmembranePotentialLR91 *mpV;
        SodiumCurrent  *mpINa;
        //SodiumCurrentLR91  *mpINa;
        BackgroundCurrentLR91 *mpIB;
        CalciumConcentrationLR91 *mpCaI;
        PlateauPotassiumCurrentLR91 *mpIKp;
        SlowInwardCurrentLR91 *mpISi;
        //PotassiumTimeInpendententCurrentLR91 *mpIK1;
        //PotassiumTimeDependentCurrentLR91 *mpIK;
       
          
>>>>>>> .r74
    public:
        // Constructor
        LR91OdeFun();
        // Destructor
        ~LR91OdeFun();
        
        // This mehtod will compute the RHS of the Luo--Rudy model
<<<<<<< .mine
        void EvaluateYDiffs(double rTime, Vec rY, Vec rYNew);   
=======
        // void EvaluateYDiffs(double rTime, Vec rY, Vec rYNew);   
        //void ComputingRHS(double tOfStimulus, std::vector<double> pInitCond); 
        std::vector<double> EvaluateYDerivatives(const double &rTime, std::vector<double> &rY);
>>>>>>> .r74
};
 
#endif //_LR91ODEFUN_HPP_
