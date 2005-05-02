#ifndef _TESTGAUSSIANQUADRATURERULE_HPP_
#define _TESTGAUSSIANQUADRATURERULE_HPP_


#include <cxxtest/TestSuite.h>
#include "GaussianQuadratureRule.hpp"



class TestGaussianQuadratureRule : public CxxTest::TestSuite 
{

public :

	/**
	 * Check points and weights are in the right ranges
	 */
	void testGaussianQuadratureRule()
	{
		// 0d
		GaussianQuadratureRule<0> quadRule01(1);
		TS_ASSERT_EQUALS(quadRule01.GetNumQuadPoints(),1);

		GaussianQuadratureRule<0> quadRule02(2);
		TS_ASSERT_EQUALS(quadRule02.GetNumQuadPoints(),1);
		
		TS_ASSERT_DELTA(quadRule01.GetWeight(0),1,1e-12);
		
		// 1d
		for (int nPoints=1; nPoints<4; nPoints++)
		{
			GaussianQuadratureRule<1> quadRule(nPoints);
			
			for(int i=0; i<quadRule.GetNumQuadPoints(); i++)
			{
				TS_ASSERT_LESS_THAN_EQUALS( quadRule.GetWeight(i),1);
				TS_ASSERT_LESS_THAN_EQUALS(-quadRule.GetWeight(i),0);
	
				TS_ASSERT_LESS_THAN( quadRule.GetQuadPoint(i)[0],1); // x<1
				TS_ASSERT_LESS_THAN(-quadRule.GetQuadPoint(i)[0],0); // x>0
			}
		}		
		
		// 2d
		for (int nPoints=1; nPoints<4; nPoints++)
		{
			GaussianQuadratureRule<2> quadRule(nPoints);
			
			for(int i=0; i<quadRule.GetNumQuadPoints(); i++)
			{
				TS_ASSERT_LESS_THAN_EQUALS( quadRule.GetWeight(i),1);
				TS_ASSERT_LESS_THAN_EQUALS(-quadRule.GetWeight(i),0);

				TS_ASSERT_LESS_THAN(-(1-quadRule.GetQuadPoint(i)[0]
									   -quadRule.GetQuadPoint(i)[1]),0); // 1-x-y>0
				TS_ASSERT_LESS_THAN(-quadRule.GetQuadPoint(i)[0],0);  // x>0
				TS_ASSERT_LESS_THAN(-quadRule.GetQuadPoint(i)[1],0);  // y>0
			}
		}

		//TODO: 3d
		
		//TODO: test through integration
	}
	
	/**
	 * Test by integrating polynomials of the form x^n.
	 * 
	 * 1d case.
	 */
	void testGaussianQuadratureRuleIntegralOneD( void )
	{
		for (int num_quad_points=1; num_quad_points<4; num_quad_points++)
		{
			GaussianQuadratureRule<1> quad_rule(num_quad_points);

			for (int poly_degree=0; poly_degree<2*num_quad_points; poly_degree++)
			{
				
				std::vector<const Node<1>*> nodes2;
				nodes2.push_back(new Node<1>(0, false, 1.0));
				nodes2.push_back(new Node<1>(1, false, 3.0));
				Element<1,1> element(nodes2);
				
				double integral=0;
				double jacobian_determinant = element.GetJacobianDeterminant();
				
		        // This assumes linear basis functions in 1d
		        double x1 = element.GetNodeLocation(0,0);
		        double x2 = element.GetNodeLocation(1,0);
		        
				for (int quad_index=0; quad_index<quad_rule.GetNumQuadPoints(); quad_index++)
				{
					Point<1> quad_point=quad_rule.GetQuadPoint(quad_index);
					Point<1> transformed_quad_point =
							Point<1>((1-quad_point[0])*x1 + quad_point[0]*x2);
							
					double integrand_value = pow(transformed_quad_point[0],poly_degree);
					
					integral+= integrand_value*jacobian_determinant
											*quad_rule.GetWeight(quad_index);
				}

				TS_ASSERT_DELTA(integral,
								1.0/(poly_degree+1.0)*(pow(3,poly_degree+1)-1),
								1e-7);
				
			}
		}
			
	}

	/**
	 * Test by integrating polynomials up to degree 2p-2, where p is the no. of
	 * points in each dimension. This is the 2d case.
	 * 
	 * We integrate things like x, y, x^2, xy, y^2, ...
	 */
	void testGaussianQuadratureRuleIntegralTwoD( void )
	{
		// Expected answers [degree_x][degree_y]
		double expected[5][5] = { {1.0/2,  1.0/6,  1.0/12, 1.0/20, 1.0/30},
								  {1.0/6,  1.0/24, 1.0/60, 1.0/120, 0},
								  {1.0/12, 1.0/60, 1.0/180, 0, 0},
								  {1.0/20, 1.0/120, 0,      0, 0},
								  {1.0/30, 0,       0,      0, 0} };
		
		for (int num_quad_points=1; num_quad_points<4; num_quad_points++)
		{
			GaussianQuadratureRule<2> quad_rule(num_quad_points);

			for (int poly_degree_x=0; poly_degree_x<2*num_quad_points-1;
				 poly_degree_x++)
			{

				for (int poly_degree_y=0;
					 poly_degree_y<2*num_quad_points-1-poly_degree_x;
					 poly_degree_y++)
				{			
					double integral = 0.0;
					
					for (int quad_index=0;
						 quad_index<quad_rule.GetNumQuadPoints();
						 quad_index++)
					{
						Point<2> quad_point=quad_rule.GetQuadPoint(quad_index);
						
						integral += pow(quad_point[0], poly_degree_x)
									*pow(quad_point[1], poly_degree_y)
									*quad_rule.GetWeight(quad_index);
					}
					
					TS_ASSERT_DELTA(integral,
									expected[poly_degree_x][poly_degree_y],
									1e-7);
				}
			}
		}
			//TS_TRACE("here gauss quad\n");
	}
		
};

#endif //_TESTGAUSSIANQUADRATURERULE_HPP_
