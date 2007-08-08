#ifndef _LINEARBASISFUNCTION_HPP_
#define _LINEARBASISFUNCTION_HPP_

#include "ChastePoint.hpp"

template <unsigned ELEM_DIM>
class LinearBasisFunction
{
public:
    static double ComputeBasisFunction(const ChastePoint<ELEM_DIM> &rPoint, unsigned basisIndex);
    static c_vector<double, ELEM_DIM> ComputeBasisFunctionDerivative(const ChastePoint<ELEM_DIM> &rPoint, unsigned basisIndex);
    
    static c_vector<double, ELEM_DIM+1> ComputeBasisFunctions(const ChastePoint<ELEM_DIM> &rPoint);
    static c_matrix<double, ELEM_DIM, ELEM_DIM+1> ComputeBasisFunctionDerivatives(const ChastePoint<ELEM_DIM> &rPoint);
    
    static c_matrix<double, ELEM_DIM, ELEM_DIM+1> ComputeTransformedBasisFunctionDerivatives(
            const ChastePoint<ELEM_DIM> &rPoint,
            const c_matrix<double, ELEM_DIM, ELEM_DIM> &rInverseJacobian);
};

/**
 * We need to specialise for the 0d case, because 0x0 matrices don't work.
 */
template <>
class LinearBasisFunction<0>
{
public:
    static double ComputeBasisFunction(const ChastePoint<0> &rPoint, unsigned basisIndex);
    static c_vector<double, 1>       ComputeBasisFunctions(const ChastePoint<0> &rPoint);
};

#endif //_LINEARBASISFUNCTION_HPP_
