/**
 * Concrete Ode4 class
 */
#include "Ode4.hpp"

Ode4::Ode4() : AbstractOdeSystem()
{
// Use AbstractOdeSystem constructors
    mNumberOfStateVariables=1;
    mInitialConditions.push_back(0.5);
}

std::vector<double> Ode4::EvaluateYDerivatives (double time, const std::vector<double> &rY)
{
    std::vector<double> y_derivatives(GetNumberOfStateVariables());
    double alpha = 100;
    y_derivatives[0]=alpha*rY[0]*(1-rY[0])*time;
    return y_derivatives;
}
