#ifndef _ABSTRACTBASISFUNCTION_HPP_
#define _ABSTRACTBASISFUNCTION_HPP_

#include "Point.hpp"
#include <vector>
#include "VectorDouble.hpp"
#include "MatrixDouble.hpp"


template <int ELEM_DIM>
class AbstractBasisFunction
{

public:
   virtual double ComputeBasisFunction(Point<ELEM_DIM> point, int basisIndex)=0;
   virtual VectorDouble ComputeBasisFunctionDerivative(Point<ELEM_DIM> point, int basisIndex)=0;
   virtual std::vector<double>       ComputeBasisFunctions(Point<ELEM_DIM> point)=0;
   virtual std::vector<VectorDouble> ComputeBasisFunctionDerivatives(Point<ELEM_DIM> point)=0;
   virtual std::vector<VectorDouble> ComputeTransformedBasisFunctionDerivatives(Point<ELEM_DIM> point, MatrixDouble inverseJacobian)=0;    
};



#endif //_ABSTRACTBASISFUNCTION_HPP_
