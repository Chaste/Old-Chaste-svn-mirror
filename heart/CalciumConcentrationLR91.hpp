#ifndef _CALCIUMCONCENTRATIONLR91_HPP_
#define _CALCIUMCONCENTRATIONLR91_HPP_


#include "IonicConcentrationLR91.hpp"
#include "ConstantsLR91.hpp"
#include <iostream>

class CalciumConcentrationLR91 : public IonicConcentrationLR91
{
    private:
        // intracellular calcium concentration 
        double mCaI;
          
    public:
        // Constructor
        CalciumConcentrationLR91();
        // Destructor
        ~CalciumConcentrationLR91();
        double GetCaI();
        void SetMagnitudeOfIonicConcentration(double caI);
        double ComputeCalciumPrime(double voltage, double d, double f, double caI, double iSi);
};


#endif //_CALCIUMCONCENTRATIONLR91_HPP_
