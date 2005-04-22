#ifndef _ABSTRACTINTEGRAND_HPP_
#define _ABSTRACTINTEGRAND_HPP_

#include "Point.hpp"

template <int SPACE_DIM>
class AbstractIntegrand
{
	public:
	virtual double Evaluate(Point<SPACE_DIM> evaluation_point)=0;
};

#endif //_ABSTRACTINTEGRAND_HPP_
