#ifndef _TESTCONSTDIRICHLETBOUNDARYCONDITION_HPP_
#define _TESTCONSTDIRICHLETBOUNDARYCONDITION_HPP_

#include <cxxtest/TestSuite.h>
#include "ConstDirichletBoundaryCondition.hpp"


class TestConstDirichletBoundaryCondition : public CxxTest::TestSuite 
{
public:
	void testConstDirichletBoundaryCondition()
	{
		Point<1> zero(0);
		ConstDirichletBoundaryCondition<1> constDirich1(1.0);
		TS_ASSERT_DELTA(constDirich1.GetValue(zero),1.0,1e-12);		
	}
};
#endif //_TESTCONSTDIRICHLETBOUNDARYCONDITION_HPP_
