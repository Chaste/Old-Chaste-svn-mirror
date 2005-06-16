#ifndef _TESTCONSTBOUNDARYCONDITION_HPP_
#define _TESTCONSTBOUNDARYCONDITION_HPP_

#include <cxxtest/TestSuite.h>
#include "ConstBoundaryCondition.hpp"
#include "VectorDouble.hpp"

class TestConstDirichletBoundaryCondition : public CxxTest::TestSuite 
{
public:
	void TestConstDirichletBoundaryConditionMethod()
	{
		Point<1> zero(0);
		ConstBoundaryCondition<1> const_dirich1(1.0);
		VectorDouble value = const_dirich1.GetValue(zero);
		TS_ASSERT_DELTA(value(0),1.0,1e-12);		
	}
};

#endif //_TESTCONSTBOUNDARYCONDITION_HPP_
