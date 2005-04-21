#ifndef _POTASSIUMTIMEDEPENDENTCURRRENTLR91_HPP_
#define _POTASSIUMTIMEDEPENDENTCURRRENTLR91_HPP_

#include "IonicCurrent.hpp"
#include "ConstantsLR91.hpp"
#include <iostream>

class PotassiumTimeDependentCurrentLR91 : public IonicCurrent
{
    private:
        // x activation gate
        double mX;
        double mAlphaX;
        double mBetaX;
        double mXi;
          
    public:
        // Constructor
        PotassiumTimeDependentCurrentLR91();
        // Destructor
        ~PotassiumTimeDependentCurrentLR91();
        double GetX();
        void UpdateAlphaAndBeta(double voltage);
        void SetGatingVariables(double x);
        void UpdateMagnitudeOfCurrent(double voltage, double x);
        double ComputeXPrime(double voltage, double x);
        void UpdateXi(double voltage);   
        double GetXi(double voltage);  
};

#endif //_POTASSIUMTIMEDEPENDENTCURRRENTLR91_HPP_

