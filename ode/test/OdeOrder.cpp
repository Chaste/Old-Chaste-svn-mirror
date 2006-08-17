// OdeOrder.cpp

#include "OdeOrder.hpp"

OdeOrder::OdeOrder() : AbstractOdeSystem()
{
// Use AbstractOdeSystem constructors
    mNumberOfStateVariables=1;
    mInitialConditions.push_back(1.0);
}

std::vector<double> OdeOrder::EvaluateYDerivatives (double time, const std::vector<double> &rY)
{
    std::vector<double> y_derivatives(GetNumberOfStateVariables());
    y_derivatives[0]=rY[0];
    return y_derivatives;
}

