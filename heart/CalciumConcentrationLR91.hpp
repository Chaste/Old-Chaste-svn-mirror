#ifndef _CALCIUMCONCENTRATIONLR91_HPP_
#define _CALCIUMCONCENTRATIONLR91_HPP_

#include "IonicConcentration.hpp"
#include "ConstantsLR91.hpp"
#include <iostream>

/**
 * Luo-Rudy 91 intracellular calcium concentration, [Ca_i].
 */
class CalciumConcentrationLR91 : public IonicConcentration
{
    public:
        // Constructor
        CalciumConcentrationLR91(void);
        // Destructor
        ~CalciumConcentrationLR91(void);
        double ComputeCalciumPrime(double voltage, double d, double f, double caI, double iSi);
};

#endif //_CALCIUMCONCENTRATIONLR91_HPP_
