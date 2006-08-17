/**
 * Concrete OdeOrderSystemOf3 class
 */
#ifndef _ODEORDERSYSTEMOF3_HPP_
#define _ODEORDERSYSTEMOF3_HPP_
#include "AbstractOdeSystem.hpp"


class OdeOrderSystemOf3 : public AbstractOdeSystem
{
public :
    OdeOrderSystemOf3();
    
    std::vector<double> EvaluateYDerivatives (double time, const std::vector<double> &rY);
    
};

#endif //_ODEORDERSYSTEMOF3_HPP_
