// TestOde3.cpp

#include "TestOde3.hpp"

TestOde3::TestOde3() : AbstractOdeSystem(2)
{
// Use AbstractOdeSystem constructors

}

void TestOde3::EvaluateYDerivatives(double t, double *y, double * yDerivatives)
{
	yDerivatives[0]=y[0]*t;
	yDerivatives[1]=y[1]*t;
}
	
