/**
 * Concrete TestOde4 class
 */ 
#include "TestOde4.hpp"

TestOde4::TestOde4() : AbstractOdeSystem(1)
{
// Use AbstractOdeSystem constructors

}

std::vector<double> TestOde4::EvaluateYDerivatives (double time, const std::vector<double> &rY)
{
	std::vector<double> yDerivatives(mNumberOfEquations);
	double alpha = 100;
	yDerivatives[0]=alpha*rY[0]*(1-rY[0])*time;
	return yDerivatives;
}
