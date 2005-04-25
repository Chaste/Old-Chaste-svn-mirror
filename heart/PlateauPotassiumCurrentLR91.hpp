#ifndef _PLATEAUPOTASSIUMCURRENTLR91_HPP_
#define _PLATEAUPOTASSIUMCURRENTLR91_HPP_

#include "IonicCurrent.hpp"
#include "ConstantsLR91.hpp"
#include <iostream>

/**
 * LR91 Plateau potassium current, IKp.
 */
class PlateauPotassiumCurrentLR91 : public IonicCurrent
{
    private:
    	// Kp is a voltage dependent term in the equation for IKp.
        double mKp;
          
    public:
        // Constructor
        PlateauPotassiumCurrentLR91();
        // Destructor
        ~PlateauPotassiumCurrentLR91();
        void UpdateMagnitudeOfCurrent(double voltage);
        void UpdateKP(double voltage);
};

#endif //_PLATEAUPOTASSIUMCURRENTLR91_HPP_
