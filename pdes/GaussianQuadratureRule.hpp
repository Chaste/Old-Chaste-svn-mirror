#ifndef _GAUSSIANQUADRATURERULE_HPP_
#define _GAUSSIANQUADRATURERULE_HPP_

#include <vector>
#include "Point.hpp"
#include "AbstractIntegrand.hpp"
#include <iostream>
#include <cmath>

template<int ELEM_DIM>
class GaussianQuadratureRule
{
	int mNumQuadPoints;	
	std::vector<double>            mWeights;
	std::vector<Point<ELEM_DIM> >  mPoints;
	
public:

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
					case 1: 
					mWeights.push_back(1);
					mPoints.push_back(Point<ELEM_DIM>(0.5)); //check
					break;
					
					case 2: 
					mWeights.push_back(0.5);
					mWeights.push_back(0.5);
					
					mPoints.push_back(Point<ELEM_DIM>(0.21132486540519));
					mPoints.push_back(Point<ELEM_DIM>(0.78867513459481));
					break;
					
					case 3: 
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
				//TODO: fill in
				assert(0);
				
				switch(numPointsInEachDimension)
				{
					case 1: 
					mWeights.push_back(1);
					mPoints.push_back(Point<ELEM_DIM>(1));
					break;
					
					case 2: 
					mWeights.push_back(1);
					mWeights.push_back(1);
					mWeights.push_back(1);
					mWeights.push_back(1);
					
					mPoints.push_back(Point<ELEM_DIM>(1));
					mPoints.push_back(Point<ELEM_DIM>(1));
					mPoints.push_back(Point<ELEM_DIM>(1));
					mPoints.push_back(Point<ELEM_DIM>(1));
					break;
					
					case 3: 
					mWeights.push_back(1);
					mWeights.push_back(1);
					mWeights.push_back(1);
					mWeights.push_back(1);
					mWeights.push_back(1);
					mWeights.push_back(1);
					mWeights.push_back(1);
					mWeights.push_back(1);
					mWeights.push_back(1);
					
					mPoints.push_back(Point<ELEM_DIM>(1));
					mPoints.push_back(Point<ELEM_DIM>(1));
					mPoints.push_back(Point<ELEM_DIM>(1));
					mPoints.push_back(Point<ELEM_DIM>(1));
					mPoints.push_back(Point<ELEM_DIM>(1));
					mPoints.push_back(Point<ELEM_DIM>(1));
					mPoints.push_back(Point<ELEM_DIM>(1));
					mPoints.push_back(Point<ELEM_DIM>(1));
					mPoints.push_back(Point<ELEM_DIM>(1));
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
	
	Point<ELEM_DIM> GetQuadPoint(int index) const
	{
		assert(index < mNumQuadPoints);
		return mPoints[index];
	}
	
	double GetWeight(int index) const
	{
		assert(index < mNumQuadPoints);
		return mWeights[index];
	}
	
	int GetNumQuadPoints() const
	{
		return mNumQuadPoints;
	}
	
	// 
	
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
