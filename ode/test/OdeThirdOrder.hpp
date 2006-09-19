/**
 * Concrete OdeThirdOrder class
 */
#ifndef _ODETHIRDORDER_HPP
#define _ODETHIRDORDER_HPP
#include "AbstractOdeSystem.hpp"


class OdeThirdOrder : public AbstractOdeSystem
{
public :
    OdeThirdOrder()
        : AbstractOdeSystem(3) // 3 here is the number of unknowns
    {
        mInitialConditions.push_back(0.0);
        mInitialConditions.push_back(1.0);
        mInitialConditions.push_back(0.0);
    }
    
    
    std::vector<double> EvaluateYDerivatives (double time, const std::vector<double> &rY)
    {
        std::vector<double> y_derivatives(GetNumberOfStateVariables());
        y_derivatives[0]=rY[0]-rY[1]+rY[2];
        y_derivatives[1]=rY[1]-rY[2];
        y_derivatives[2]=2*rY[1]-rY[2];
        return y_derivatives;
    }
};

#endif //_ODETHIRDORDER_HPP
