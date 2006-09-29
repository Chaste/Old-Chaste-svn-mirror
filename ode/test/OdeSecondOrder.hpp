#ifndef _ODESECONDORDER_HPP_
#define _ODESECONDORDER_HPP_
#include "AbstractOdeSystem.hpp"


class OdeSecondOrder : public AbstractOdeSystem
{
public :
    OdeSecondOrder()
            : AbstractOdeSystem(2) // 2 here is the number of unknowns
    {
        mInitialConditions.push_back(0.0);
        mInitialConditions.push_back(1.0);
    }
    
    
    std::vector<double> EvaluateYDerivatives (double time, const std::vector<double> &rY)
    {
        std::vector<double> y_derivatives(GetNumberOfStateVariables());
        y_derivatives[0] =  rY[1];
        y_derivatives[1] = -rY[0];
        return y_derivatives;
    }
};

#endif //_ODESECONDORDER_HPP_
