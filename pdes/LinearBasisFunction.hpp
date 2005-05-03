#ifndef _LINEARBASISFUNCTION_HPP_
#define _LINEARBASISFUNCTION_HPP_

#include "AbstractBasisFunction.hpp"
#include "Point.hpp"

template <int ELEM_DIM>
class LinearBasisFunction : public AbstractBasisFunction<ELEM_DIM>
{
    
public:

    double ComputeBasisFunction(Point<ELEM_DIM> point, int basisIndex) const;
    VectorDouble ComputeBasisFunctionDerivative(Point<ELEM_DIM> point, int basisIndex) const;
    
    std::vector<double>       ComputeBasisFunctions(Point<ELEM_DIM> psi) const;
    std::vector<VectorDouble> ComputeBasisFunctionDerivatives(Point<ELEM_DIM> psi) const;
 
 	std::vector<VectorDouble> ComputeTransformedBasisFunctionDerivatives(Point<ELEM_DIM> psi, 
 	                                                                     MatrixDouble inverseJacobian) const;
    
};

#endif //_LINEARBASISFUNCTION_HPP_
