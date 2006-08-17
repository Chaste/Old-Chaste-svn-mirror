#ifndef _ODEORDER_HPP_
#define _ODEORDER_HPP_
#include "AbstractOdeSystem.hpp"


class OdeOrder : public AbstractOdeSystem
{
public :
    OdeOrder();
    
    std::vector<double> EvaluateYDerivatives (double time, const std::vector<double> &rY);
    
};

#endif //_ODEORDER_HPP_
