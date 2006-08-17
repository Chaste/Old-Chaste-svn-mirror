/**
 * Concrete Ode4 class
 */
#ifndef _ODE4_HPP_
#define _ODE4_HPP_
#include "AbstractOdeSystem.hpp"


class Ode4 : public AbstractOdeSystem
{
public :
    Ode4();
    
    std::vector<double> EvaluateYDerivatives (double time, const std::vector<double> &rY);
    
};

#endif //_ODE4_HPP_

