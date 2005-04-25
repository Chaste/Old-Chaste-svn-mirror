#ifndef _TESTCONSTBOUNDARYCONDITION_HPP_
#define _TESTCONSTBOUNDARYCONDITION_HPP_

#include <cxxtest/TestSuite.h>
#include "ConstBoundaryCondition.hpp"


class TestConstDirichletBoundaryCondition : public CxxTest::TestSuite 
{
public:
	void testConstDirichletBoundaryCondition()
	{
		Point<1> zero(0);
		ConstBoundaryCondition<1> constDirich1(1.0);
		TS_ASSERT_DELTA(constDirich1.GetValue(zero),1.0,1e-12);		
	}
};
#endif //_TESTCONSTBOUNDARYCONDITION_HPP_
