#ifndef _ABSTRACTBOUNDARYCONDITION_HPP_
#define _ABSTRACTBOUNDARYCONDITION_HPP_

#include "Point.hpp"

/**
 * Abstract base class for boundary conditions.
 */
template<int SPACE_DIM>
class AbstractBoundaryCondition
{
public:
	/**
	 * Get the value of the boundary condition at a given point.
	 * 
	 * @param x The point at which to evaluate the boundary condition.
	 */
	virtual double GetValue(const Point<SPACE_DIM> x) const = 0;
};

#endif //_ABSTRACTDIRICHLETBOUNDARYCONDITION_HPP_
