// TestOde1.cpp

#include "TestOde1.hpp"

TestOde1::TestOde1() : AbstractOdeSystem(1)
{
// Use AbstractOdeSystem constructors

}

void TestOde1::EvaluateYDerivatives(double t, double *y, double * yDerivatives)
{
	yDerivatives[0]=t;
}
	
