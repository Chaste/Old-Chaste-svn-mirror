#ifndef _TESTODE2_HPP_
#define _TESTODE2_HPP_
// TestOde2.hpp
#include "AbstractOdeSystem.hpp"


class TestOde2 : public AbstractOdeSystem
{
	public :
	TestOde2();
	
	void EvaluateYDerivatives (double rTime, double * rY, double * rYDerivatives);

};

#endif //_TESTODE2_HPP_
