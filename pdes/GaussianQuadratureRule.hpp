#ifndef _GAUSSIANQUADRATURERULE_HPP_
#define _GAUSSIANQUADRATURERULE_HPP_

#include <vector>
#include "Point.hpp"

template<int ELEM_DIM>
class GaussianQuadratureRule
{
	int mNumQuadPoints;	
	std::vector<double>            mWeights;
	std::vector<Point<ELEM_DIM> >  mPoints;
	
public:

	GaussianQuadratureRule(unsigned int numPointsInEachDimension)
	{
		assert(numPointsInEachDimension >  0);
		assert(numPointsInEachDimension <= 3);
		
		mNumQuadPoints = numPointsInEachDimension^(ELEM_DIM);	

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
	
	int GetNumQuadPoints()
	{
		return mNumQuadPoints;
	}
	
};





#endif //_GAUSSIANQUADRATURERULE_HPP_
