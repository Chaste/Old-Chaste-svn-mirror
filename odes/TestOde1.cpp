// TestOde1.cpp

#include "TestOde1.hpp"

TestOde1::TestOde1() : AbstractOdeSystem(1)
{
// Use AbstractOdeSystem constructors

}

std::vector<double> TestOde1::EvaluateYDerivatives (double rTime, std::vector<double> &rY)
{
	std::vector<double> rYDerivatives(mNumberOfEquations);
	rYDerivatives[0]=1.0;
	
	return rYDerivatives;
}
	
