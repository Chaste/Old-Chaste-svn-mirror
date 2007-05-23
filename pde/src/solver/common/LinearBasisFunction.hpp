#ifndef _LINEARBASISFUNCTION_HPP_
#define _LINEARBASISFUNCTION_HPP_

#include "Point.hpp"

template <unsigned ELEM_DIM>
class LinearBasisFunction
{
public:
    static double ComputeBasisFunction(const Point<ELEM_DIM> &rPoint, unsigned basisIndex);
    static c_vector<double, ELEM_DIM> ComputeBasisFunctionDerivative(const Point<ELEM_DIM> &rPoint, unsigned basisIndex);
    
    static c_vector<double, ELEM_DIM+1> ComputeBasisFunctions(const Point<ELEM_DIM> &rPoint);
    static c_matrix<double, ELEM_DIM, ELEM_DIM+1> ComputeBasisFunctionDerivatives(const Point<ELEM_DIM> &rPoint);
    
    static c_matrix<double, ELEM_DIM, ELEM_DIM+1> ComputeTransformedBasisFunctionDerivatives(
            const Point<ELEM_DIM> &rPoint,
            const c_matrix<double, ELEM_DIM, ELEM_DIM> &rInverseJacobian);
};

/**
 * We need to specialise for the 0d case, because 0x0 matrices don't work.
 */
template <>
class LinearBasisFunction<0>
{
public:
    static double ComputeBasisFunction(const Point<0> &rPoint, unsigned basisIndex);
    static c_vector<double, 1>       ComputeBasisFunctions(const Point<0> &rPoint);
};

#endif //_LINEARBASISFUNCTION_HPP_
