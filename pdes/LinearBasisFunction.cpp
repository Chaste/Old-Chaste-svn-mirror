#include "LinearBasisFunction.hpp"
#include "Point.hpp"
#include <cassert>

template <int ELEM_DIM>
double LinearBasisFunction<ELEM_DIM>::ComputeBasisFunction(Point<ELEM_DIM> point, int basisIndex)
{
	switch(ELEM_DIM)
	{
	case 1:
		assert(basisIndex == 0 || basisIndex == 1);
    	if(basisIndex == 0)
    	{	
        	return 1.0 - point[0];
    	}
    	else
    	{
        	return point[0];   
    	}
    	break;
    	
    case 2:
    	assert(basisIndex == 0 || basisIndex == 1 || basisIndex == 2);
    	if(basisIndex == 0)
    	{	
        	return 1.0 - point[0] - point[1];
    	}
    	else if(basisIndex == 1)
    	{	
        	return point[0];
    	}
    	else
    	{
        	return point[1];   
    	}
    	break;
    	
    	case 3:
    	assert(basisIndex == 0 || basisIndex == 1 || basisIndex == 2 || basisIndex == 3);
    	if(basisIndex == 0)
    	{	
        	return 1.0 - point[0] - point[1] - point[2];
    	}
    	else if(basisIndex == 1)
    	{	
        	return point[0];
    	}
    	else if(basisIndex == 2)
    	{	
        	return point[1];
    	}
    	else
    	{
        	return point[2];   
    	}
    	break;
	}
}


template <int ELEM_DIM>
VectorDouble LinearBasisFunction<ELEM_DIM>::ComputeBasisFunctionDerivative(Point<ELEM_DIM> point, int basisIndex)
{
    VectorDouble gradN(ELEM_DIM);
    
    if(basisIndex == 0)
    {
        gradN(0) = -1;
    }
    else
    {
        gradN(0) = 1;
    }
    
    return gradN;
}



template <int ELEM_DIM>
std::vector<double> LinearBasisFunction<ELEM_DIM>::ComputeBasisFunctions(Point<ELEM_DIM> psi)
{
    std::vector<double> basisValues(ELEM_DIM+1);
    
    for(int i=0;i<ELEM_DIM+1;i++)
    {
        basisValues[i] = ComputeBasisFunction(psi,i);
    }
    return basisValues;
}



template <int ELEM_DIM>
std::vector<VectorDouble>  LinearBasisFunction<ELEM_DIM>::ComputeBasisFunctionDerivatives(Point<ELEM_DIM> psi)
{
    std::vector<VectorDouble> basisGradValues;
    
    for(int i=0;i<ELEM_DIM+1;i++)
    {
        basisGradValues.push_back(ComputeBasisFunctionDerivative(psi,i));
    }

    return basisGradValues;    
}
