/*
 * Concrete Ode2 class
 */
#ifndef _ODE2_HPP_
#define _ODE2_HPP_
#include "AbstractOdeSystem.hpp"


class Ode2 : public AbstractOdeSystem
{
public :

    Ode2() : AbstractOdeSystem(1) // 1 here is the number of unknowns
    {
        mInitialConditions.push_back(4.0);
    }
    
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
    {
        rDY[0]=rY[0]*time;
    }
};

#endif //_ODE2_HPP_
