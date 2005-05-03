#ifndef _ABSTRACTNONLINEARELLIPTICPDE_HPP_
#define _ABSTRACTNONLINEARELLIPTICPDE_HPP_

class MatrixDouble;
#include "Point.hpp"

template <int SPACE_DIM>
class AbstractNonlinearEllipticPde
{
public:

    virtual double ComputeLinearSourceTerm(Point<SPACE_DIM> x)=0;
    
    virtual double ComputeNonlinearSourceTerm(Point<SPACE_DIM> x,
                                              double u)=0;

    virtual MatrixDouble ComputeDiffusionTerm(Point<SPACE_DIM> x,
                                              double u)=0;
                                              
	virtual MatrixDouble ComputeDiffusionTermPrime(Point<SPACE_DIM> x,
                                              double u)=0;
                                              
    virtual double ComputeNonlinearSourceTermPrime(Point<SPACE_DIM> x,
                                              double u)=0;                                          

};

#endif //_ABSTRACTNONLINEARELLIPTICPDE_HPP_
