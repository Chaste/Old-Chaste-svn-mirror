#ifndef _ABSTRACTBASISFUNCTION_HPP_
#define _ABSTRACTBASISFUNCTION_HPP_

#include "Point.hpp"
#include <vector>
#include "VectorDouble.hpp"
#include "MatrixDouble.hpp"

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
   virtual void                      ComputeBasisFunctionsWithUpdate(const Point<ELEM_DIM> &rPoint, std::vector<double> &rBasisValues) const =0;
   virtual std::vector<VectorDouble> ComputeBasisFunctionDerivatives(const Point<ELEM_DIM> &rPoint) const =0;
   virtual std::vector<c_vector<double, ELEM_DIM> > ComputeTransformedBasisFunctionDerivatives(const Point<ELEM_DIM> &rPoint, const c_matrix<double, ELEM_DIM, ELEM_DIM> &rInverseJacobian) const =0;
   virtual ~AbstractBasisFunction() { };
};

/**
 * We need to specialise for the 0d case, because 0x0 matrices don't work.
 */
template <>
class AbstractBasisFunction<0>
{
public:
   virtual double ComputeBasisFunction(const Point<0> &rPoint, int basisIndex) const =0;
   virtual std::vector<double> ComputeBasisFunctions(const Point<0> &rPoint) const =0;
   virtual void ComputeBasisFunctionsWithUpdate(const Point<0> &rPoint, std::vector<double> &rBasisValues) const =0;
   virtual ~AbstractBasisFunction() { };
};

#endif //_ABSTRACTBASISFUNCTION_HPP_
