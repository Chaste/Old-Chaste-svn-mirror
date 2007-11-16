#ifndef _ABSTRACTBOUNDARYCONDITION_HPP_
#define _ABSTRACTBOUNDARYCONDITION_HPP_

#include "ChastePoint.hpp"
#include "UblasCustomFunctions.hpp"

/**
 * Abstract base class for boundary conditions.
 */
template<unsigned SPACE_DIM>
class AbstractBoundaryCondition
{
public:
    /**
     * Get the value of the boundary condition at a given point.
     * 
     * @param x The point at which to evaluate the boundary condition.
     */
    virtual double GetValue(const ChastePoint<SPACE_DIM>& x) const = 0;
    
    // Make derived classes work
    AbstractBoundaryCondition()
    {}
    virtual ~AbstractBoundaryCondition()
    {}
};

#endif //_ABSTRACTDIRICHLETBOUNDARYCONDITION_HPP_
