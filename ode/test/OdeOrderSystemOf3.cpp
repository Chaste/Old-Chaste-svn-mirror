/**
 * Concrete OdeOrderSystemOf3 class
 */ 
#include "OdeOrderSystemOf3.hpp"

OdeOrderSystemOf3::OdeOrderSystemOf3() : AbstractOdeSystem()
{
// Use AbstractOdeSystem constructors
    mNumberOfStateVariables=3;
    mInitialConditions.push_back(0.0);
    mInitialConditions.push_back(1.0);
    mInitialConditions.push_back(0.0);
}

std::vector<double> OdeOrderSystemOf3::EvaluateYDerivatives (double time, const std::vector<double> &rY)
{
	std::vector<double> y_derivatives(GetNumberOfStateVariables());
	y_derivatives[0]=rY[0]-rY[1]+rY[2];
	y_derivatives[1]=rY[1]-rY[2];
	y_derivatives[2]=2*rY[1]-rY[2];
	return y_derivatives;
}
	
