#ifndef _PLATEAUPOTASSIUMCURRENTLR91_HPP_
#define _PLATEAUPOTASSIUMCURRENTLR91_HPP_

#include "IonicCurrent.hpp"
#include "ConstantsLR91.hpp"
#include <iostream>

class PlateauPotassiumCurrentLR91 : public IonicCurrent
{
    private:
        // gating variable
        double mKp;
          
    public:
        // Constructor
        PlateauPotassiumCurrentLR91();
        // Destructor
        ~PlateauPotassiumCurrentLR91();
        double GetKp();
        void UpdateMagnitudeOfCurrent(double voltage);
        void UpdateGatingVariables(double voltage);
};

#endif //_PLATEAUPOTASSIUMCURRENTLR91_HPP_
