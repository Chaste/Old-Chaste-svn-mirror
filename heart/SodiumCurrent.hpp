#ifndef _SODIUMCURRENT_HPP_
#define _SODIUMCURRENT_HPP_

#include "IonicCurrent.hpp"
#include "ConstantsLR91.hpp"
#include <iostream>

class SodiumCurrent : public IonicCurrent
{
    private:
        // m activation gate
        double mM;
        // h,j inactivation gates
        double mH;
        double mJ;
        double mAlphaM;
        double mAlphaH;
        double mAlphaJ;
        double mBetaM;
        double mBetaH;
        double mBetaJ;
          
    public:
        // Constructor
        SodiumCurrent();
        // Destructor
        ~SodiumCurrent();
        double GetM();
        double GetH();
        double GetJ();
        void UpdateMagnitudeOfCurrent(double voltage, double m, double h, double j);
        void UpdateAlphaAndBeta(double voltage);
        void UpdateGatingVariables(double m, double h, double j);
        double ComputeMPrime(double voltage, double m, double h, double j);
        double ComputeHPrime(double voltage, double m, double h, double j);
        double ComputeJPrime(double voltage, double m, double h, double j);
};

#endif //_SODIUMCURRENT_HPP_
