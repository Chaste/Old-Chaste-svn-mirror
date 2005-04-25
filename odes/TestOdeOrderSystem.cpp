// TestOdeOrderSystem.cpp

#include "TestOdeOrderSystem.hpp"
#include <math.h>

TestOdeOrderSystem::TestOdeOrderSystem() : AbstractOdeSystem(2)
{
// Use AbstractOdeSystem constructors

}

std::vector<double> TestOdeOrderSystem::EvaluateYDerivatives (double time, const std::vector<double> &rY)
{
	std::vector<double> yDerivatives(mNumberOfEquations);
	yDerivatives[0]=rY[1];
	yDerivatives[1]=-rY[0];
	return yDerivatives;
}
	
