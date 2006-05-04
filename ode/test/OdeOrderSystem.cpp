// OdeOrderSystem.cpp

#include "OdeOrderSystem.hpp"

OdeOrderSystem::OdeOrderSystem() : AbstractOdeSystem()
{
// Use AbstractOdeSystem constructors
    mNumberOfStateVariables=2;
    mInitialConditions.push_back(0.0);
    mInitialConditions.push_back(1.0);
}

std::vector<double> OdeOrderSystem::EvaluateYDerivatives (double time, const std::vector<double> &rY)
{
	std::vector<double> y_derivatives(GetNumberOfStateVariables());
	y_derivatives[0]=rY[1];
	y_derivatives[1]=-rY[0];
	return y_derivatives;
}
	
