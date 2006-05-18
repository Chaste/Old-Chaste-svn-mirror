#ifndef _LINEARBASISFUNCTION_HPP_
#define _LINEARBASISFUNCTION_HPP_

#include "AbstractBasisFunction.hpp"
#include "Point.hpp"

template <int ELEM_DIM>
class LinearBasisFunction : public AbstractBasisFunction<ELEM_DIM>
{
public:
    double ComputeBasisFunction(const Point<ELEM_DIM> &rPoint, int basisIndex) const;
    c_vector<double, ELEM_DIM> ComputeBasisFunctionDerivative(const Point<ELEM_DIM> &rPoint, int basisIndex) const;
    
    std::vector<double>       ComputeBasisFunctions(const Point<ELEM_DIM> &rPoint) const;
    void                      ComputeBasisFunctionsWithUpdate(const Point<ELEM_DIM> &rPoint, std::vector<double> &rBasisValues) const;
    c_matrix<double, ELEM_DIM, ELEM_DIM+1> ComputeBasisFunctionDerivatives(const Point<ELEM_DIM> &rPoint) const;
 
 	c_matrix<double, ELEM_DIM, ELEM_DIM+1> ComputeTransformedBasisFunctionDerivatives(const Point<ELEM_DIM> &rPoint,
 	                                                                     const c_matrix<double, ELEM_DIM, ELEM_DIM> &rInverseJacobian) const;
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
