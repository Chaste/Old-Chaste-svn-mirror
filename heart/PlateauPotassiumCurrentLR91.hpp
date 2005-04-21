#ifndef _PLATEAUPOTASSIUMCURRENTLR91_HPP_
#define _PLATEAUPOTASSIUMCURRENTLR91_HPP_

#include "IonicCurrent.hpp"
#include "ConstantsLR91.hpp"
#include <iostream>

class PlateauPotassiumCurrentLR91 : public IonicCurrent
{
    private:
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
