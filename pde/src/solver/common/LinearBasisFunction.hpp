#ifndef _LINEARBASISFUNCTION_HPP_
#define _LINEARBASISFUNCTION_HPP_

#include "AbstractBasisFunction.hpp"
#include "Point.hpp"

template <int ELEM_DIM>
class LinearBasisFunction : public AbstractBasisFunction<ELEM_DIM>
{
    
public:

    double ComputeBasisFunction(const Point<ELEM_DIM> &rPoint, int basisIndex) const;
    VectorDouble ComputeBasisFunctionDerivative(const Point<ELEM_DIM> &rPoint, int basisIndex) const;
    
    std::vector<double>       ComputeBasisFunctions(const Point<ELEM_DIM> &rPoint) const;
    std::vector<VectorDouble> ComputeBasisFunctionDerivatives(const Point<ELEM_DIM> &rPoint) const;
 
 	std::vector<VectorDouble> ComputeTransformedBasisFunctionDerivatives(const Point<ELEM_DIM> &rPoint,
 	                                                                     const MatrixDouble &rInverseJacobian) const;
    
};

#endif //_LINEARBASISFUNCTION_HPP_
