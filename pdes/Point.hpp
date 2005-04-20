#ifndef _POINT_HPP_
#define _POINT_HPP_

#include <cassert>

template<int DIM> 
class Point
{
private:
    double mLocation[DIM];

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

    double operator[] (unsigned i) const
    {
        assert(i<DIM); 
        return mLocation[i];
    }
    
    void SetCoordinate(unsigned i, double value)
    {
         assert(i<DIM);
         mLocation[i] = value;   
    } 
    
};



#endif //_POINT_HPP_
