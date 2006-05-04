/**
 * Concrete Ode1 class
 */ 
#include "Ode1.hpp"

Ode1::Ode1() : AbstractOdeSystem()
{
// Use AbstractOdeSystem constructors
    mNumberOfStateVariables=1;
    mInitialConditions.push_back(0.0);
}

std::vector<double> Ode1::EvaluateYDerivatives (double time, const std::vector<double> &rY)
{
	std::vector<double> y_derivatives(GetNumberOfStateVariables());
	y_derivatives[0]=1.0;
	
	return y_derivatives;
}
	
