// TestOde2.cpp

#include "TestOde2.hpp"

TestOde2::TestOde2() : AbstractOdeSystem(1)
{
// Use AbstractOdeSystem constructors

}

void TestOde2::EvaluateYDerivatives(double t, double *y, double * yDerivatives)
{
	yDerivatives[0]=y[0]*t;
}
	
