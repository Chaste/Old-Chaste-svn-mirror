#ifndef _TESTODE1_HPP_
#define _TESTODE1_HPP_
// TestOde1.hpp
#include "AbstractOdeSystem.hpp"


class TestOde1 : public AbstractOdeSystem
{
	public :
	TestOde1();
	
	void EvaluateYDerivatives (double rTime, double * rY, double * rYDerivatives);

};

#endif //_TESTODE1_HPP_
