#ifndef _QUADRATICBASISFUNCTION_HPP_
#define _QUADRATICBASISFUNCTION_HPP_

#include "AbstractBasisFunction.hpp"
#include "Point.hpp"
template <int ELEM_DIM>

class QuadraticBasisFunction : public AbstractBasisFunction<ELEM_DIM>
{
	public:

    double ComputeBasisFunction(const Point<ELEM_DIM> &rPoint, int basisIndex) const;
    VectorDouble ComputeBasisFunctionDerivative(const Point<ELEM_DIM> &rPoint, int basisIndex) const;
    
    std::vector<double>       ComputeBasisFunctions(const Point<ELEM_DIM> &rPoint) const;
    void                      ComputeBasisFunctionsWithUpdate(const Point<ELEM_DIM> &rPoint, std::vector<double> &rBasisValues) const;
    std::vector<VectorDouble> ComputeBasisFunctionDerivatives(const Point<ELEM_DIM> &rPoint) const;
 
 	std::vector<VectorDouble> ComputeTransformedBasisFunctionDerivatives(const Point<ELEM_DIM> &rPoint,
 	                                                                     const MatrixDouble &rInverseJacobian) const;
    
	

};

#endif //_QUADRATICBASISFUNCTION_HPP_
