#ifndef _TESTODEORDERSYSTEM_HPP_
#define _TESTODEORDERSYSTEM_HPP_
// TestOdeOrderSystem.hpp
#include "AbstractOdeSystem.hpp"


class TestOdeOrderSystem : public AbstractOdeSystem
{
	public :
	TestOdeOrderSystem();
	
	std::vector<double> TestOdeOrderSystem::EvaluateYDerivatives (double rTime, std::vector<double> &rY);

};

#endif //_TESTODEORDERSYSTEM_HPP_
