#ifndef _GAUSSIANQUADRATURERULE_HPP_
#define _GAUSSIANQUADRATURERULE_HPP_

#include <vector>
#include "Point.hpp"
#include "AbstractIntegrand.hpp"
#include <iostream>
#include <cmath>

/**
 * This class encapsulates tables of gaussian quadrature points and the
 * associated weights.
 * 
 * Data is available for 1d and 2d quadrature over (canonical) triangles,
 * with between 1 and 3 (inclusive) gauss points in each dimension.
 * 
 * TODO: Add 3d.
 */

template<int ELEM_DIM>
class GaussianQuadratureRule
{
	int mNumQuadPoints;	
	std::vector<double>            mWeights;
	std::vector<Point<ELEM_DIM> >  mPoints;
	
public:

	/**
	 * The constructor builds the appropriate table for the dimension (given
	 * by the template argument) and number of points in each dimension (given
	 * as a constructor argument).
	 */
	GaussianQuadratureRule(int numPointsInEachDimension)
	{
		assert(numPointsInEachDimension >  0);
		assert(numPointsInEachDimension <= 3);
		
		mNumQuadPoints = (int) pow((double) numPointsInEachDimension,(ELEM_DIM));	

		mWeights.reserve(mNumQuadPoints);
		mPoints.reserve(mNumQuadPoints);

		switch(ELEM_DIM)
		{
			case 1 :
			{
				switch(numPointsInEachDimension)
				{
					case 1: // 1d, 1 point
					mWeights.push_back(1);
					mPoints.push_back(Point<ELEM_DIM>(0.5)); //check
					break;
					
					case 2: // 1d, 2 points
					mWeights.push_back(0.5);
					mWeights.push_back(0.5);
					
					mPoints.push_back(Point<ELEM_DIM>(0.21132486540519));
					mPoints.push_back(Point<ELEM_DIM>(0.78867513459481));
					break;
					
					case 3: // 1d, 3 points
					mWeights.push_back(5.0/18.0);
					mWeights.push_back(4.0/9.0);
					mWeights.push_back(5.0/18.0);
					
					mPoints.push_back(Point<ELEM_DIM>(0.1127016654));
					mPoints.push_back(Point<ELEM_DIM>(0.5));
					mPoints.push_back(Point<ELEM_DIM>(0.8872983346));
					break;
				}
			}
			break;
			case 2 :
			{				
				switch(numPointsInEachDimension)
				{
					case 1: // 2d, 1 point per dimension
					mWeights.push_back(0.5);
					mPoints.push_back(Point<ELEM_DIM>(0.25,0.5));
					break;
					
					case 2: // 2d, 2 points per dimension
					mWeights.push_back(0.19716878364870);
					mWeights.push_back(0.19716878364870);
					mWeights.push_back(0.05283121635130);
					mWeights.push_back(0.05283121635130);
					
					mPoints.push_back(Point<ELEM_DIM>(0.16666666666667,0.21132486540519));
					mPoints.push_back(Point<ELEM_DIM>(0.62200846792815,0.21132486540519));
					mPoints.push_back(Point<ELEM_DIM>(0.04465819873852,0.78867513459481));
					mPoints.push_back(Point<ELEM_DIM>(0.16666666666667,0.78867513459481));
					break;
					
					case 3: // 2d, 3 points per dimension
					mWeights.push_back(0.06846437766975);
					mWeights.push_back(0.10954300427160);
					mWeights.push_back(0.06846437766975);
					mWeights.push_back(0.06172839506173);
					mWeights.push_back(0.09876543209877);
					mWeights.push_back(0.06172839506173);
					mWeights.push_back(0.00869611615741);
					mWeights.push_back(0.01391378585185);
					mWeights.push_back(0.00869611615741);
					
					mPoints.push_back(Point<ELEM_DIM>(0.10000000001607,0.11270166540000));
					mPoints.push_back(Point<ELEM_DIM>(0.44364916730000,0.11270166540000));
					mPoints.push_back(Point<ELEM_DIM>(0.78729833458393,0.11270166540000));
					mPoints.push_back(Point<ELEM_DIM>(0.05635083270000,0.50000000000000));
					mPoints.push_back(Point<ELEM_DIM>(0.25000000000000,0.50000000000000));
					mPoints.push_back(Point<ELEM_DIM>(0.44364916730000,0.50000000000000));
					mPoints.push_back(Point<ELEM_DIM>(0.01270166538393,0.88729833460000));
					mPoints.push_back(Point<ELEM_DIM>(0.05635083270000,0.88729833460000));
					mPoints.push_back(Point<ELEM_DIM>(0.10000000001607,0.88729833460000));
					break;
				}
				
			}
			break;
			case 3 :
			{
				//TODO: fill in
				assert(0);				
			}
			break;
			
			default:
			// TODO: something more sensible
			assert(0);
		}

	}
	
	/**
	 * Get a quadrature point.
	 * 
	 * @param index The index of the point to return.
	 * @return A gaussian quadrature point.
	 */
	Point<ELEM_DIM> GetQuadPoint(int index) const
	{
		assert(index < mNumQuadPoints);
		return mPoints[index];
	}
	
	/**
	 * Get the weight associated with a quadrature point.
	 */
	double GetWeight(int index) const
	{
		assert(index < mNumQuadPoints);
		return mWeights[index];
	}
	
	/**
	 * Get the number of quadrature points. This is the number of points in 
	 * each dimension, raised to the power of the number of dimensions.
	 */
	int GetNumQuadPoints() const
	{
		return mNumQuadPoints;
	}
	
	
	/**
	 * This method is only used by TestGaussianQuadratureRule at the moment.
	 * 
	 * We assume ELEM_DIM=SPACE_DIM
	 * TODO: Integrate when ELEM_DIM<SPACE_DIM ?
	 */
	double Integrate(Element<ELEM_DIM,ELEM_DIM> &rElement,
						AbstractIntegrand<ELEM_DIM> &rFunction,
						AbstractIntegrand<ELEM_DIM> &rCanonicalFunction) const
	{
		double integral=0;
		
		double jacobian_determinant = rElement.GetJacobianDeterminant();
		
        // This assumes linear basis functions in 1d
        double x1 = rElement.GetNodeLocation(0,0);
        double x2 = rElement.GetNodeLocation(1,0);
        
 
		
		for (int quad_index=0; quad_index<GetNumQuadPoints(); quad_index++)
		{
			Point<ELEM_DIM> quad_point=GetQuadPoint(quad_index);
			// TODO: extend above 1d
			Point<ELEM_DIM> transformed_quad_point =
					Point<ELEM_DIM>((1-quad_point[0])*x1 + quad_point[0]*x2);
			double integrand_value=
					rFunction.Evaluate(transformed_quad_point)
					* rCanonicalFunction.Evaluate(quad_point);
								
					integral+= integrand_value*jacobian_determinant
									*GetWeight(quad_index);
		}
		return integral;
	}
	
};





#endif //_GAUSSIANQUADRATURERULE_HPP_
