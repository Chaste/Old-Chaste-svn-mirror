// TestOdeOrder.cpp

#include "TestOdeOrder.hpp"
#include <math.h>

TestOdeOrder::TestOdeOrder() : AbstractOdeSystem(1)
{
// Use AbstractOdeSystem constructors

}

std::vector<double> TestOdeOrder::EvaluateYDerivatives (double time, const std::vector<double> &rY)
{
	std::vector<double> yDerivatives(mNumberOfEquations);
	yDerivatives[0]=rY[0];
	return yDerivatives;
}
	
