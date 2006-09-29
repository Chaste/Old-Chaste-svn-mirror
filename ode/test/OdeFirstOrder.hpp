#ifndef _ODEFIRSTORDER_HPP_
#define _ODEFIRSTORDER_HPP_
#include "AbstractOdeSystem.hpp"


class OdeFirstOrder : public AbstractOdeSystem
{
public :
    OdeFirstOrder()
            : AbstractOdeSystem(1) // 1 here is the number of variables
    {
        mInitialConditions.push_back(1.0);
    }
    
    
    std::vector<double> EvaluateYDerivatives (double time, const std::vector<double> &rY)
    {
        std::vector<double> y_derivatives(GetNumberOfStateVariables());
        y_derivatives[0] = rY[0];
        return y_derivatives;
    }
    
};

#endif //_ODEFIRSTORDER_HPP_
