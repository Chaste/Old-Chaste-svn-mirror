/**
 * Concrete TestOde3 class
 */ 
#include "TestOde3.hpp"

TestOde3::TestOde3() : AbstractOdeSystem()
{
// Use AbstractOdeSystem constructors
    mNumberOfStateVariables=2;
    mInitialConditions.push_back(4.0);
    mInitialConditions.push_back(8.0);
}

std::vector<double> TestOde3::EvaluateYDerivatives (double time, const std::vector<double> &rY)
{
	std::vector<double> yDerivatives(GetNumberOfStateVariables());
	yDerivatives[0]=rY[0]*time;
	yDerivatives[1]=rY[1]*time;
	return yDerivatives;
}
	
