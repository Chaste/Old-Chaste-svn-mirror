/**
 * Concrete TestOde3 class
 */ 
#ifndef _TESTODE3_HPP_
#define _TESTODE3_HPP_
// TestOde3.hpp
#include "AbstractOdeSystem.hpp"


class TestOde3 : public AbstractOdeSystem
{
	public :
	TestOde3();
	
	std::vector<double> TestOde3::EvaluateYDerivatives (double rTime, std::vector<double> &rY);

};

#endif //_TESTODE3_HPP_
