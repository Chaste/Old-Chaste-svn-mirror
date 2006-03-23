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
    VectorDouble *mpValue;
    
public:
	/**
	 * Create a new boundary condition object.
	 * 
	 * @param value The value of this boundary condition at all points where it
	 *    is applied.
	 */
    ConstBoundaryCondition(const double value)
    {
    	mpValue = new VectorDouble(1);
    	(*mpValue)(0) = value;
    }

    ConstBoundaryCondition(const VectorDouble value)
    {
    	mpValue = new VectorDouble( value.Size() );
    	for( int i = 0; i < value.Size(); i++ )
    	{
    		(*mpValue)(i) = value(i);
    	}
    }
    
    /**
     * @param x The point at which this boundary condition is to be evaluated.
     * @return The constant value given in the constructor.
     */
    VectorDouble GetValue( const Point<SPACE_DIM> ) const
    {
        return *mpValue;
    }
    
    /**
     * Delete the stored VectorDouble giving the value of this boundary condition.
     */
    ~ConstBoundaryCondition()
    {
    	delete mpValue;
    }
};

#endif //_CONSTBOUNDARYCONDITION_HPP_
