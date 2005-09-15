// TestOdeOrderSystem.cpp

#include "TestOdeOrderSystem.hpp"
#include <math.h>

TestOdeOrderSystem::TestOdeOrderSystem() : AbstractOdeSystem()
{
// Use AbstractOdeSystem constructors

    mInitialConditions.push_back(0.0);
    mInitialConditions.push_back(1.0);
}

std::vector<double> TestOdeOrderSystem::EvaluateYDerivatives (double time, const std::vector<double> &rY)
{
	std::vector<double> yDerivatives(GetNumberOfStateVariables());
	yDerivatives[0]=rY[1];
	yDerivatives[1]=-rY[0];
	return yDerivatives;
}
	
