/**
 * Concrete TestOdeOrderSystemOf3 class
 */ 
#include "TestOdeOrderSystemOf3.hpp"
#include <math.h>

TestOdeOrderSystemOf3::TestOdeOrderSystemOf3() : AbstractOdeSystem()
{
// Use AbstractOdeSystem constructors
    mNumberOfStateVariables=3;
    mInitialConditions.push_back(0.0);
    mInitialConditions.push_back(1.0);
    mInitialConditions.push_back(0.0);
}

std::vector<double> TestOdeOrderSystemOf3::EvaluateYDerivatives (double time, const std::vector<double> &rY)
{
	std::vector<double> yDerivatives(GetNumberOfStateVariables());
	yDerivatives[0]=rY[0]-rY[1]+rY[2];
	yDerivatives[1]=rY[1]-rY[2];
	yDerivatives[2]=2*rY[1]-rY[2];
	return yDerivatives;
}
	
