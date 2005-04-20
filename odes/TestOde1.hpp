#ifndef _TESTODE1_HPP_
#define _TESTODE1_HPP_
// TestOde1.hpp
#include "AbstractOdeSystem.hpp"


class TestOde1 : public AbstractOdeSystem
{
	public :
	TestOde1(double t, double * y);
	
	void EvaluateYPrime (double rTime, double * rY, double * rYPrime);

};

#endif //_TESTODE1_HPP_
