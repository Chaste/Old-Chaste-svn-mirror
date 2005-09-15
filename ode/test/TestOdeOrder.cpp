// TestOdeOrder.cpp

#include "TestOdeOrder.hpp"
#include <math.h>

TestOdeOrder::TestOdeOrder() : AbstractOdeSystem()
{
// Use AbstractOdeSystem constructors

    mInitialConditions.push_back(1.0);
}

std::vector<double> TestOdeOrder::EvaluateYDerivatives (double time, const std::vector<double> &rY)
{
	std::vector<double> yDerivatives(GetNumberOfStateVariables());
	yDerivatives[0]=rY[0];
	return yDerivatives;
}
	
