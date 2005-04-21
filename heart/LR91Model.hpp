#ifndef _LR91MODEL_HPP_
#define _LR91MODEL_HPP_

#include "LR91OdeFun.hpp"
#include "AbstractStimulusFunction.hpp"
#include <iostream>

class LR91Model
{
    private:
        LR91OdeFun *mpLR91OdeSystem;
        double mV;
        double mM;
        double mH;
        double mJ;
        double mD;
        double mF;
        double mX;
        double mCaI;
    
        
    public:
        // Constructor
        // initializes the LR91Model with initial conditions!
        LR91Model(double voltage, double m, double h, double j, double d, 
                  double f, double x, double caI, AbstractStimulusFunction *pStimulus);
        // Destructor
        ~LR91Model();
        
        //Solve should solve LR91 system and return whatever type the solver returns
        std::vector<double> Solve();//OdeSolver myOdeSolver);
        
       
};
 

#endif //_LR91MODEL_HPP_
