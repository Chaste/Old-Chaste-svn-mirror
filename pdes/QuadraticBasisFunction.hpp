#ifndef _QUADRATICBASISFUNCTION_HPP_
#define _QUADRATICBASISFUNCTION_HPP_

#include "AbstractBasisFunction.hpp"
#include "Point.hpp"
template <int ELEM_DIM>

class QuadraticBasisFunction : public AbstractBasisFunction<ELEM_DIM>
{
	public:

    double ComputeBasisFunction(Point<ELEM_DIM> point, int basisIndex);
    VectorDouble ComputeBasisFunctionDerivative(Point<ELEM_DIM> point, int basisIndex);
    
    std::vector<double>       ComputeBasisFunctions(Point<ELEM_DIM> psi);
    std::vector<VectorDouble> ComputeBasisFunctionDerivatives(Point<ELEM_DIM> psi);
 
 	std::vector<VectorDouble> ComputeTransformedBasisFunctionDerivatives(Point<ELEM_DIM> psi, 
 	                                                                     MatrixDouble inverseJacobian);
    
	

};

#endif //_QUADRATICBASISFUNCTION_HPP_
