#ifndef _LINEARBASISFUNCTION_HPP_
#define _LINEARBASISFUNCTION_HPP_

#include "AbstractBasisFunction.hpp"
#include "Point.hpp"

template <int ELEM_DIM>
class LinearBasisFunction : public AbstractBasisFunction<ELEM_DIM>
{
    
public:

    double ComputeBasisFunction(Point<ELEM_DIM> point, int basisIndex);
    VectorDouble ComputeBasisFunctionDerivative(Point<ELEM_DIM> point, int basisIndex);
    
    std::vector<double>       ComputeBasisFunctions(Point<ELEM_DIM> psi);
    std::vector<VectorDouble> ComputeBasisFunctionDerivatives(Point<ELEM_DIM> psi);
 
    
};



#endif //_LINEARBASISFUNCTION_HPP_
