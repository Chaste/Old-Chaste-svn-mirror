#ifndef _SLOWINWARDCURRENTLR91_HPP_
#define _SLOWINWARDCURRENTLR91_HPP_

#include "IonicCurrent.hpp"
#include "ConstantsLR91.hpp"
#include <iostream>
/**
 *LR91 Slow Inward Current ISi.
 */
class SlowInwardCurrentLR91 : public IonicCurrent
{
    private:
        // gating variables
        double mD;
        double mF;
        // Nernst Potential for Calcium 
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
        //Update Nerst potential 
        void UpdateESi(double caI);
        double ComputeDPrime(double voltage, double d, double f);
        double ComputeFPrime(double voltage, double d, double f);        
};

#endif //_SLOWINWARDCURRENTLR91_HPP_
