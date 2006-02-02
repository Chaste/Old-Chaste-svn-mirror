#ifndef _ABSTRACTBASISFUNCTION_HPP_
#define _ABSTRACTBASISFUNCTION_HPP_

#include "Point.hpp"
#include <vector>
class VectorDouble;
class MatrixDouble;

/**
 * Abstract base class for basis functions. There are methods to compute
 * the value and derivative of a particular basis function, or all basis
 * functions on an element together.
 * 
 * The methods are documented more fully in the LinearBasisFunction class.
 * 
 * @see LinearBasisFunction
 */
template <int ELEM_DIM>
class AbstractBasisFunction
{

public:
   virtual double ComputeBasisFunction(const Point<ELEM_DIM> &rPoint, int basisIndex) const =0;
   virtual VectorDouble ComputeBasisFunctionDerivative(const Point<ELEM_DIM> &rPoint, int basisIndex) const =0;
   virtual std::vector<double>       ComputeBasisFunctions(const Point<ELEM_DIM> &rPoint) const =0;
   virtual std::vector<VectorDouble> ComputeBasisFunctionDerivatives(const Point<ELEM_DIM> &rPoint) const =0;
   virtual std::vector<VectorDouble> ComputeTransformedBasisFunctionDerivatives(const Point<ELEM_DIM> &rPoint, const MatrixDouble &rInverseJacobian) const =0;
};



#endif //_ABSTRACTBASISFUNCTION_HPP_
