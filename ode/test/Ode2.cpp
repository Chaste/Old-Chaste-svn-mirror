/**
 * Concrete Ode2 class
 */
#include "Ode2.hpp"

Ode2::Ode2() : AbstractOdeSystem()
{
// Use AbstractOdeSystem constructors
    mNumberOfStateVariables=1;
    mInitialConditions.push_back(4.0);
}

std::vector<double> Ode2::EvaluateYDerivatives (double time, const std::vector<double> &rY)
{
    std::vector<double> y_derivatives(GetNumberOfStateVariables());
    y_derivatives[0]=rY[0]*time;
    return y_derivatives;
}

