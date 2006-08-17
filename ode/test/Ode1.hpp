/**
 * Concrete Ode1 class
 */
#ifndef _ODE1_HPP_
#define _ODE1_HPP_
#include "AbstractOdeSystem.hpp"


class Ode1 : public AbstractOdeSystem
{
public :
    Ode1();
    
    std::vector<double> EvaluateYDerivatives (double time, const std::vector<double> &rY);
    
};

#endif //_ODE1_HPP_
