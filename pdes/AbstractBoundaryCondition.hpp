#ifndef _ABSTRACTBOUNDARYCONDITION_HPP_
#define _ABSTRACTBOUNDARYCONDITION_HPP_

#include "Point.hpp"

template<int SPACE_DIM>
class AbstractBoundaryCondition
{
public:
    virtual double GetValue(Point<SPACE_DIM> x) = 0;
};


#endif //_ABSTRACTDIRICHLETBOUNDARYCONDITION_HPP_
