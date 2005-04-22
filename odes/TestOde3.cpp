/**
 * Concrete TestOde3 class
 */ 
#include "TestOde3.hpp"

TestOde3::TestOde3() : AbstractOdeSystem(2)
{
// Use AbstractOdeSystem constructors

}

std::vector<double> TestOde3::EvaluateYDerivatives (double rTime, std::vector<double> &rY)
{
	std::vector<double> yDerivatives(mNumberOfEquations);
	yDerivatives[0]=rY[0]*rTime;
	yDerivatives[1]=rY[1]*rTime;
	return yDerivatives;
}
	
