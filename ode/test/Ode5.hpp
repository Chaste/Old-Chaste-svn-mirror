#ifndef ODE5_HPP_
#define ODE5_HPP_

#include "AbstractOdeSystem.hpp"


class Ode5 : public AbstractOdeSystem
{
public :
    Ode5()
            : AbstractOdeSystem(1)  // 1 here is the number of unknowns
    {
        mInitialConditions.push_back(0.2);
    }
    
    std::vector<double> EvaluateYDerivatives (double time, const std::vector<double> &rY)
    {
        std::vector<double> y_derivatives(GetNumberOfStateVariables());
        double alpha = 100;
        y_derivatives[0]=alpha*rY[0]*(1-rY[0]);
        return y_derivatives;
    }
    
};

#endif /*ODE5_HPP_*/
