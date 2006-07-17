#ifndef _POINT_HPP_
#define _POINT_HPP_


#include <boost/numeric/ublas/vector.hpp>
#include <cassert>
#include <vector>
using namespace boost::numeric::ublas;

template<int DIM> 
class Point
{
private:
  c_vector<double, DIM> mLocation;

public:

  /**
   * Create a Point object.
   * There are 3 optional arguments, which can be used to specify the values
   * of the first 3 dimensions, if present.
   * 
   * Point now uses a ublas vector to store its location. The
   * rGetLocation method returns a reference to this vector.
   * Use of this method together with ublas operations
   * is the perfered way to use this class.
   */
  Point(double v1=0, double v2=0, double v3=0)
  {
    if (DIM>0)
    {
      mLocation[0] = v1;
    }
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
      mLocation(i) = coords.at(i);
    }
  }
    
  c_vector<double, DIM>& rGetLocation(void)
  {
    return mLocation;
  }

  double operator[] (unsigned i) const
  {
    assert((int)i<DIM); 
    return mLocation(i);
  }
    
  void SetCoordinate(unsigned i, double value)
  {
    assert(i<DIM);
    mLocation(i) = value;   
  }
  
  Point<DIM> MidPoint(Point<DIM> otherPoint)
  {
    Point<DIM> mid_point;
    c_vector<double, DIM>& mid_point_location = mid_point.rGetLocation();
    mid_point_location = 0.5 * (mLocation + otherPoint.rGetLocation());
    return mid_point;
  }     
};

template<>
class Point<0>
{
public:

  /**
   * Create a zero-dimensional Point object.
   * There are 3 optional arguments, which should not be used.
   */
  Point(double v1=0, double v2=0, double v3=0)
  {
  }
};


#endif //_POINT_HPP_
