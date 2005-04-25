#ifndef _CONSTDIRICHLETBOUNDARYCONDITION_HPP_
#define _CONSTDIRICHLETBOUNDARYCONDITION_HPP_

#include "AbstractBoundaryCondition.hpp"

template<int SPACE_DIM>
class ConstDirichletBoundaryCondition : public AbstractBoundaryCondition<SPACE_DIM>
{	
private :
	double mValue;
	
public :
	ConstDirichletBoundaryCondition(double value) : mValue(value) {};
	
	double GetValue( const Point<SPACE_DIM> x) const
    {
   		return mValue;
    }
};

#endif //_CONSTDIRICHLETBOUNDARYCONDITION_HPP_
