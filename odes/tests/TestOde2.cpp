/**
 * Concrete TestOde2 class
 */
#include "TestOde2.hpp"

TestOde2::TestOde2() : AbstractOdeSystem(1)
{
// Use AbstractOdeSystem constructors

}

std::vector<double> TestOde2::EvaluateYDerivatives (double time, const std::vector<double> &rY)
{
	std::vector<double> yDerivatives(mNumberOfEquations);
	yDerivatives[0]=rY[0]*time;
	return yDerivatives;
}
	
