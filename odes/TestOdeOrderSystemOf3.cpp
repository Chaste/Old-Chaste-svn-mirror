/**
 * Concrete TestOdeOrderSystemOf3 class
 */ 
#include "TestOdeOrderSystemOf3.hpp"
#include <math.h>

TestOdeOrderSystemOf3::TestOdeOrderSystemOf3() : AbstractOdeSystem(3)
{
// Use AbstractOdeSystem constructors

}

std::vector<double> TestOdeOrderSystemOf3::EvaluateYDerivatives (double rTime, std::vector<double> &rY)
{
	std::vector<double> yDerivatives(mNumberOfEquations);
	yDerivatives[0]=rY[0]-rY[1]+rY[2];
	yDerivatives[1]=rY[1]-rY[2];
	yDerivatives[2]=2*rY[1]-rY[2];
	return yDerivatives;
}
	
