/**
 * Concrete TestOde1 class
 */ 
#include "TestOde1.hpp"

TestOde1::TestOde1() : AbstractOdeSystem(1)
{
// Use AbstractOdeSystem constructors

}

std::vector<double> TestOde1::EvaluateYDerivatives (double time, const std::vector<double> &rY)
{
	std::vector<double> yDerivatives(GetNumberOfEquations());
	yDerivatives[0]=1.0;
	
	return yDerivatives;
}
	
