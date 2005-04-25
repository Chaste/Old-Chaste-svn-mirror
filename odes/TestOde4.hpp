/**
 * Concrete TestOde4 class
 */ 
#ifndef _TESTODE4_HPP_
#define _TESTODE4_HPP_
// TestOde4.hpp
#include "AbstractOdeSystem.hpp"


class TestOde4 : public AbstractOdeSystem
{
	public :
	TestOde4();
	
	std::vector<double> TestOde4::EvaluateYDerivatives (double time, const std::vector<double> &rY);

};

#endif //_TESTODE4_HPP_

