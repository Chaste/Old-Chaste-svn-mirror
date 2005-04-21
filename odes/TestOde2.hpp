#ifndef _TESTODE2_HPP_
#define _TESTODE2_HPP_
// TestOde2.hpp
#include "AbstractOdeSystem.hpp"


class TestOde2 : public AbstractOdeSystem
{
	public :
	TestOde2();
	
	std::vector<double> TestOde2::EvaluateYDerivatives (double rTime, std::vector<double> &rY);

};

#endif //_TESTODE2_HPP_
