/**
 * Concrete TestOdeOrderSystemOf3 class
 */ 
#ifndef _TESTODEORDERSYSTEMOF3_HPP_
#define _TESTODEORDERSYSTEMOF3_HPP_
// TestOdeOrderSystemOf3.hpp
#include "AbstractOdeSystem.hpp"


class TestOdeOrderSystemOf3 : public AbstractOdeSystem
{
	public :
	TestOdeOrderSystemOf3();
	
	std::vector<double> TestOdeOrderSystemOf3::EvaluateYDerivatives (double rTime, std::vector<double> &rY);

};

#endif //_TESTODEORDERSYSTEMOF3_HPP_
