#ifndef _CONSTBOUNDARYCONDITION_HPP_
#define _CONSTBOUNDARYCONDITION_HPP_

#include "AbstractBoundaryCondition.hpp"

/**
 * Boundary condition that takes a constant value wherever it is applied.
 */
template<int SPACE_DIM>
class ConstBoundaryCondition : public AbstractBoundaryCondition<SPACE_DIM>
{   
private:
    double mValue;
    
public:
	/**
	 * Create a new boundary condition object.
	 * 
	 * @param value The value of this boundary condition at all points where it
	 *    is applied.
	 */
    ConstBoundaryCondition(double value) : mValue(value) {};
    
    /**
     * @param x The point at which this boundary condition is to be evaluated.
     * @return The constant value given in the constructor.
     */
    double GetValue( const Point<SPACE_DIM> x) const
    {
        return mValue;
    }
};

#endif //_CONSTBOUNDARYCONDITION_HPP_
