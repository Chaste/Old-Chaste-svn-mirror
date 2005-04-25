#ifndef _TESTODEORDER_HPP_
#define _TESTODEORDER_HPP_
// TestOdeOrder.hpp
#include "AbstractOdeSystem.hpp"


class TestOdeOrder : public AbstractOdeSystem
{
	public :
	TestOdeOrder();
	
	std::vector<double> TestOdeOrder::EvaluateYDerivatives (double time, const std::vector<double> &rY);

};

#endif //_TESTODEORDER_HPP_
