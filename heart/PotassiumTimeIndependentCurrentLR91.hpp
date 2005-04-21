#ifndef _POTASSIUMTIMEINDEPENDENTCURRENTLR91_HPP_
#define _POTASSIUMTIMEINDEPENDENTCURRENTLR91_HPP_
#include "IonicCurrent.hpp"
#include "ConstantsLR91.hpp"
#include <iostream>

class PotassiumTimeIndependentCurrentLR91 : public IonicCurrent
{
    private:
        double mK1; // K1 inactivation gate, do not have to solve K1'!
        double mAlphaK1;
        double mBetaK1;
          
    public:
        // Constructor
        PotassiumTimeIndependentCurrentLR91();
        // Destructor
        ~PotassiumTimeIndependentCurrentLR91();
        void UpdateAlphaAndBeta(double voltage);
        void UpdateK1(double voltage);
        void UpdateMagnitudeOfCurrent(double voltage);
};

#endif //_POTASSIUMTIMEINDEPENDENTCURRENTLR91_HPP_
