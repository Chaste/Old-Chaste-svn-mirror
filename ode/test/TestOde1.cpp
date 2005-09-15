/**
 * Concrete TestOde1 class
 */ 
#include "TestOde1.hpp"

TestOde1::TestOde1() : AbstractOdeSystem()
{
// Use AbstractOdeSystem constructors

    mInitialConditions.push_back(0.0);
}

std::vector<double> TestOde1::EvaluateYDerivatives (double time, const std::vector<double> &rY)
{
	std::vector<double> yDerivatives(GetNumberOfStateVariables());
	yDerivatives[0]=1.0;
	
	return yDerivatives;
}
	
