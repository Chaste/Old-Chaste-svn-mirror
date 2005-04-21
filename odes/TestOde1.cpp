// TestOde1.cpp

#include "TestOde1.hpp"

TestOde1::TestOde1(double t, double * y) : AbstractOdeSystem(1, t, y)
{
// Use AbstractOdeSystem constructors	
}

void TestOde1::EvaluateYDerivatives(double t, double *y, double * yDerivatives)
{
	yDerivatives[0]=t;
}
	
