#ifndef _POINT_HPP_
#define _POINT_HPP_

#include <cassert>

template<int DIM> 
class Point
{
private:
    double mLocation[DIM];

public:

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
