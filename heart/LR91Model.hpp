#ifndef _LR91MODEL_HPP_
#define _LR91MODEL_HPP_

#include "LR91OdeFun.hpp"
#include <iostream>

#include "petscvec.h"

class LR91Model
{
    private:
        LR91OdeFun myL91OdeSystem;
       // OdeSolver myOdeSolver;
        
    public:
        // Constructor
        //LR91Model(double voltage, double m, double h, double j, double d, double f, double x, double caI);
        // Destructor
        //~LR91Model();
        
        
       
};
 

#endif //_LR91MODEL_HPP_
