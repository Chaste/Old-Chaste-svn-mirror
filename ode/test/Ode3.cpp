/**
 * Concrete Ode3 class
 */
#include "Ode3.hpp"

Ode3::Ode3() : AbstractOdeSystem()
{
// Use AbstractOdeSystem constructors
    mNumberOfStateVariables=2;
    mInitialConditions.push_back(4.0);
    mInitialConditions.push_back(8.0);
}

std::vector<double> Ode3::EvaluateYDerivatives (double time, const std::vector<double> &rY)
{
    std::vector<double> y_derivatives(GetNumberOfStateVariables());
    y_derivatives[0]=rY[0]*time;
    y_derivatives[1]=rY[1]*time;
    return y_derivatives;
}

