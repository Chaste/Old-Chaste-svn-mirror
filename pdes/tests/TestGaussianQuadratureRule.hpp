#ifndef _TESTGAUSSIANQUADRATURERULE_HPP_
#define _TESTGAUSSIANQUADRATURERULE_HPP_


#include <cxxtest/TestSuite.h>
#include "GaussianQuadratureRule.hpp"
#include "PolyFunction.hpp"



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
	
	void testGaussianQuadratureRuleIntegralOneD( void )
	{
		for (int num_quad_points=1; num_quad_points<4; num_quad_points++)
		{
			GaussianQuadratureRule<1> quad_rule(num_quad_points);

			for (int poly_degree=0; poly_degree<1; poly_degree++)
			{

				poly_function function(poly_degree);
				poly_function unit_function(0);
				
				std::vector<Node<1>*> nodes;
				nodes.push_back(new Node<1>(0, false, 0.0));
				nodes.push_back(new Node<1>(1, false, 1.0));
				Element<1,1> element(nodes);				
				
				TS_ASSERT_DELTA(quad_rule.Integrate(element, function, unit_function),
								1.0/(poly_degree+1.0),
								1e-7);
				
				std::vector<Node<1>*> nodes2;
				nodes2.push_back(new Node<1>(0, false, 1.0));
				nodes2.push_back(new Node<1>(1, false, 3.0));
				Element<1,1> element2(nodes2);
				
				TS_ASSERT_DELTA(quad_rule.Integrate(element2, function, unit_function),
								1.0/(poly_degree+1.0)*(pow(3,poly_degree+1)-1),
								1e-7);
				
			}
		}
			
	}
		
};

#endif //_TESTGAUSSIANQUADRATURERULE_HPP_
