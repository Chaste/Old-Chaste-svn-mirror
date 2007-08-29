/*
 * Concrete Ode1 class
 */
#ifndef _ODE1_HPP_
#define _ODE1_HPP_
#include "AbstractOdeSystem.hpp"


class Ode1 : public AbstractOdeSystem
{
public :

    Ode1()  : AbstractOdeSystem(1) // 1 here is the number of variables
    {
        mInitialConditions.push_back(0.0);
        mStateVariables = mInitialConditions;
    }
    
    
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
    {
        rDY[0]=1.0;
    }
};

#endif //_ODE1_HPP_
