#ifndef _POTASSIUMTIMEINDEPENDENTCURRENT_HPP_
#define _POTASSIUMIMEINDEPENDENTCURRENT_HPP_
#include "IonicCurrent.hpp"
#include "ConstantsLR91.hpp"
#include <iostream>

class PotassiumTimeIndependentCurrent : public IonicCurrent
{
    private:
        // m activation gate
        double mM;
        // h,j inactivation gates
        double mH;
        double mJ;
          
    public:
        // Constructor
        PotassiumTimeIndependentCurrent(double m, double h, double j);
        // Destructor
        ~PotassiumTimeIndependentCurrent();
        double GetM();
        double GetH();
        double GetJ();
        void UpdateMagnitudeOfCurrent(double voltage);
};
#endif //_POTASSIUMTIMEINDEPENDENTCURRENT_HPP_
