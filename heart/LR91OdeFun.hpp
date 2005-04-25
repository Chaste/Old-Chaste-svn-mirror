#ifndef _LR91ODEFUN_HPP_
#define _LR91ODEFUN_HPP_

#include "ConstantsLR91.hpp"
#include "AbstractOdeSystem.hpp"
#include "TransmembranePotentialLR91.hpp"
#include "SodiumCurrentLR91.hpp"
#include "BackgroundCurrentLR91.hpp"
#include "CalciumConcentrationLR91.hpp"
#include "PlateauPotassiumCurrentLR91.hpp"
#include "SlowInwardCurrentLR91.hpp"
#include "PotassiumTimeIndependentCurrentLR91.hpp"
#include "PotassiumTimeDependentCurrentLR91.hpp"
#include "AbstractStimulusFunction.hpp"
#include <iostream>
#include <vector>

/**
 * This class sets up the Lou-Rudy 91 system of equations.
 */
class LR91OdeFun : public AbstractOdeSystem
{
    private:
    	// Current and voltage components (objects) of the LR91 model
        TransmembranePotentialLR91 *mpV;
        SodiumCurrentLR91  *mpINa;
        BackgroundCurrentLR91 *mpIB;
        CalciumConcentrationLR91 *mpCaI;
        PlateauPotassiumCurrentLR91 *mpIKp;
        SlowInwardCurrentLR91 *mpISi;
        PotassiumTimeIndependentCurrentLR91 *mpIK1;
        PotassiumTimeDependentCurrentLR91 *mpIK;
        AbstractStimulusFunction *mpStimulus;

    public:
        // Constructor
        LR91OdeFun(AbstractStimulusFunction *stimulus);
        // Destructor
        ~LR91OdeFun(void);
        
        // This method will compute the RHS of the Luo-Rudy 91 model
        std::vector<double> EvaluateYDerivatives(double time, const std::vector<double> &rY);
};
 
#endif //_LR91ODEFUN_HPP_
