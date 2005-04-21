#ifndef _TESTGAUSSIANQUADRATURERULE_HPP_
#define _TESTGAUSSIANQUADRATURERULE_HPP_


#include <cxxtest/TestSuite.h>
#include "GaussianQuadratureRule.hpp"

class TestGaussianQuadratureRule : public CxxTest::TestSuite 
{

public :
	void testGaussianQuadratureRule()
	{
		GaussianQuadratureRule<1> quadRule(3);
		
		for(int i=0; i<quadRule.GetNumQuadPoints(); i++)
		{
			
			TS_ASSERT_LESS_THAN_EQUALS( quadRule.GetWeight(i),1);
			TS_ASSERT_LESS_THAN_EQUALS(-quadRule.GetWeight(i),0);

			TS_ASSERT_LESS_THAN( quadRule.GetQuadPoint(i)[0],1);
			TS_ASSERT_LESS_THAN(-quadRule.GetQuadPoint(i)[0],0);
		}
		
				
		//TODO: 2d and 3d
		
		//TODO: test through integration
	}
};

#endif //_TESTGAUSSIANQUADRATURERULE_HPP_
