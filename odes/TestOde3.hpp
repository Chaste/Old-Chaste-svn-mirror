#ifndef _TESTODE3_HPP_
#define _TESTODE3_HPP_
// TestOde3.hpp
#include "AbstractOdeSystem.hpp"


class TestOde3 : public AbstractOdeSystem
{
	public :
	TestOde3();
	
	void EvaluateYDerivatives (double rTime, double * rY, double * rYDerivatives);

};

#endif //_TESTODE3_HPP_
