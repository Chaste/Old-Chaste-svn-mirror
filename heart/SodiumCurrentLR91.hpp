#ifndef _SODIUMCURRENTLR91_HPP_
#define _SODIUMCURRENTLR91_HPP_

#include "IonicCurrent.hpp"
#include "ConstantsLR91.hpp"
#include <iostream>

/**
 * LR91 Fast Sodium Current INa
 * 
 */

class SodiumCurrentLR91 : public IonicCurrent
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
        SodiumCurrentLR91();
        // Destructor
        ~SodiumCurrentLR91();
        double GetM();
        double GetH();
        double GetJ();
        void UpdateMagnitudeOfCurrent(double voltage, double m, double h, double j);
        void UpdateAlphaAndBeta(double voltage);
        void SetGatingVariables(double m, double h, double j);
        double ComputeMPrime(double voltage, double m, double h, double j);
        double ComputeHPrime(double voltage, double m, double h, double j);
        double ComputeJPrime(double voltage, double m, double h, double j);
};

#endif //_SODIUMCURRENTLR91_HPP_
