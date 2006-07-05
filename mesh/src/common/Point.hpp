#ifndef _POINT_HPP_
#define _POINT_HPP_

#include <cassert>
#include <vector>

template<int DIM> 
class Point
{
private:
    double mLocation[DIM?DIM:1]; // Bodge to allow zero dimensional points

public:

	/**
	 * Create a Point object.
	 * There are 3 optional arguments, which can be used to specify the values
	 * of the first 3 dimensions, if present.
	 */
	Point(double v1=0, double v2=0, double v3=0)
	{
		mLocation[0] = v1;
		if (DIM>1)
		{
			mLocation[1] = v2;
		}
		if (DIM>2)
		{
			mLocation[2] = v3;
		}		
	}
	
	/**
	 * Create a Point object.
	 * This constructor takes a vector giving the coordinates of the point.
	 * The length of the vector must be at least the dimension of the point.
	 */
	Point(std::vector<double> coords)
	{
		for (int i=0; i<DIM; i++)
		{
			mLocation[i] = coords.at(i);
		}
	}

    double operator[] (unsigned i) const
    {
        assert((int)i<DIM); 
        return mLocation[i];
    }
    
    void SetCoordinate(unsigned i, double value)
    {
         assert(i<DIM);
         mLocation[i] = value;   
    }
    
    Point<DIM> MidPoint(Point<DIM> otherPoint)
    {
    	Point<DIM> new_point;
    	for (int i=0; i<DIM; i++)
    	{
    		new_point.SetCoordinate(i, 0.5*(this->mLocation[i] + otherPoint.mLocation[i]));
    	}
    	return new_point;
    } 
    
};

#endif //_POINT_HPP_
