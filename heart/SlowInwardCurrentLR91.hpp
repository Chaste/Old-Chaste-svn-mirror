#ifndef _SLOWINWARDCURRENTLR91_HPP_
#define _SLOWINWARDCURRENTLR91_HPP_

#include "IonicCurrent.hpp"
#include "ConstantsLR91.hpp"
#include <iostream>

class SlowInwardCurrentLR91 : public IonicCurrent
{
    private:
        // gating variable
        double mD;
        double mF;
        double mESi;
        double mAlphaD;
        double mAlphaF;
        double mBetaD;
        double mBetaF;
          
    public:
        // Constructor
        SlowInwardCurrentLR91();
        // Destructor
        ~SlowInwardCurrentLR91();
        double GetD();
        double GetF();
        void UpdateMagnitudeOfCurrent(double voltage, double d, double f, double caI);
        void UpdateAlphaAndBeta(double voltage);
        void SetGatingVariables(double d, double f);
        // Update Nerst potential that depends on Cai
        void UpdateESi(double caI);
        double ComputeDPrime(double voltage, double d, double f);
        double ComputeFPrime(double voltage, double d, double f);        
};

#endif //_SLOWINWARDCURRENTLR91_HPP_
