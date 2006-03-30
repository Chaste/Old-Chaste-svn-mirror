/**
 * Concrete TestOde2 class
 */
#include "TestOde2.hpp"

TestOde2::TestOde2() : AbstractOdeSystem()
{
// Use AbstractOdeSystem constructors
    mNumberOfStateVariables=1;
    mInitialConditions.push_back(4.0);
}

std::vector<double> TestOde2::EvaluateYDerivatives (double time, const std::vector<double> &rY)
{
	std::vector<double> yDerivatives(GetNumberOfStateVariables());
	yDerivatives[0]=rY[0]*time;
	return yDerivatives;
}
	
