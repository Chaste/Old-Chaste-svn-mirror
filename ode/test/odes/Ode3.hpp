/*
 * Concrete Ode3 class
 */
#ifndef _ODE3_HPP_
#define _ODE3_HPP_
#include "AbstractOdeSystem.hpp"


class Ode3 : public AbstractOdeSystem
{
public :

    Ode3()  : AbstractOdeSystem(2) // 2 here is the number of unknowns
    {
        mInitialConditions.push_back(4.0);
        mInitialConditions.push_back(8.0);
    }
    
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
    {
        rDY[0]=rY[0]*time;
        rDY[1]=rY[1]*time;
    }
};

#endif //_ODE3_HPP_
