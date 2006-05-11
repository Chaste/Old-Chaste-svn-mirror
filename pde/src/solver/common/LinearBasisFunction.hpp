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
    void                      ComputeBasisFunctionsWithUpdate(const Point<ELEM_DIM> &rPoint, std::vector<double> &rBasisValues) const;
    std::vector<VectorDouble> ComputeBasisFunctionDerivatives(const Point<ELEM_DIM> &rPoint) const;
 
 	std::vector<c_vector<double, ELEM_DIM> > ComputeTransformedBasisFunctionDerivatives(const Point<ELEM_DIM> &rPoint,
 	                                                                     const MatrixDouble &rInverseJacobian) const;
};

/**
 * We need to specialise for the 0d case, because 0x0 matrices don't work.
 */
template <>
class LinearBasisFunction<0> : public AbstractBasisFunction<0>
{
public:
    double ComputeBasisFunction(const Point<0> &rPoint, int basisIndex) const;
    std::vector<double>       ComputeBasisFunctions(const Point<0> &rPoint) const;
    void                      ComputeBasisFunctionsWithUpdate(const Point<0> &rPoint, std::vector<double> &rBasisValues) const;
};

#endif //_LINEARBASISFUNCTION_HPP_
