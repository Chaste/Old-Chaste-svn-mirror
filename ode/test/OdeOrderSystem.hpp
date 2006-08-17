#ifndef _ODEORDERSYSTEM_HPP_
#define _ODEORDERSYSTEM_HPP_
#include "AbstractOdeSystem.hpp"


class OdeOrderSystem : public AbstractOdeSystem
{
public :
    OdeOrderSystem();
    
    std::vector<double> EvaluateYDerivatives (double time, const std::vector<double> &rY);
    
};

#endif //_ODEORDERSYSTEM_HPP_
