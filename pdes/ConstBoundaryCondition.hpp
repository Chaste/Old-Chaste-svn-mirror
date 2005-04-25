#ifndef _CONSTBOUNDARYCONDITION_HPP_
#define _CONSTBOUNDARYCONDITION_HPP_

#include "AbstractBoundaryCondition.hpp"

template<int SPACE_DIM>
class ConstBoundaryCondition : public AbstractBoundaryCondition<SPACE_DIM>
{   
private :
    double mValue;
    
public :
    ConstBoundaryCondition(double value) : mValue(value) {};
    
    double GetValue( const Point<SPACE_DIM> x) const
    {
        return mValue;
    }
};
#endif //_CONSTBOUNDARYCONDITION_HPP_
