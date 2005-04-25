#ifndef _ABSTRACTBOUNDARYCONDITION_HPP_
#define _ABSTRACTBOUNDARYCONDITION_HPP_

#include "Point.hpp"

template<int SPACE_DIM>
class AbstractBoundaryCondition
{
public:
//TODO:?	virtual ~AbstractBoundaryCondition() {}
	virtual double GetValue(const Point<SPACE_DIM> x) const = 0;
};

#endif //_ABSTRACTDIRICHLETBOUNDARYCONDITION_HPP_
