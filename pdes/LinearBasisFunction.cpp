#include "LinearBasisFunction.hpp"
#include "Point.hpp"
#include <cassert>

/**
 * Linear basis functions are implemented entirely for 1D at present.
 * 2D and 3D derivatives still need to be implemented.
 */

template <int ELEM_DIM>
double LinearBasisFunction<ELEM_DIM>::ComputeBasisFunction(Point<ELEM_DIM> point, int basisIndex)
{
	assert(ELEM_DIM < 4 && ELEM_DIM > 0);
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
    assert(ELEM_DIM < 4 && ELEM_DIM > 0);
    
	switch(ELEM_DIM)
	{
	case 1:
		assert(basisIndex == 0 || basisIndex == 1);
    	if(basisIndex == 0)
    	{	
        	gradN(0) =  - 1;
    	}
    	else
    	{
        	gradN(0) = 1;   
    	}
    	break;
    	
    case 2:
    	assert(basisIndex == 0 || basisIndex == 1 || basisIndex == 2);
    	if(basisIndex == 0)
    	{	
        	gradN(0) =  - 1;
        	gradN(1) = - 1;
    	}
    	else if(basisIndex == 1)
    	{	
        	gradN(0) = 1;
        	gradN(1) = 0;
    	}
    	else
    	{
        	gradN(0) = 0;
        	gradN(1) = 1;   
    	}
    	break;
    	
    case 3:
    	assert(basisIndex == 0 || basisIndex == 1 || basisIndex == 2 || basisIndex == 3);
    	if(basisIndex == 0)
    	{	
        	gradN(0) =  - 1;
        	gradN(1) = - 1;
        	gradN(2) = - 1;
    	}
    	else if(basisIndex == 1)
    	{	
        	gradN(0) =  1;
        	gradN(1) =  0;
        	gradN(2) =  0;
    	}
    	else if(basisIndex == 2)
    	{	
        	gradN(0) =  0;
        	gradN(1) =  1;
        	gradN(2) =  0;
    	}
    	else
    	{
        	gradN(0) =  0;
        	gradN(1) =  0;
        	gradN(2) =  1;   
    	}
    	break;
	}    
    return gradN;
}



template <int ELEM_DIM>
std::vector<double> LinearBasisFunction<ELEM_DIM>::ComputeBasisFunctions(Point<ELEM_DIM> psi)
{
    std::vector<double> basisValues(ELEM_DIM+1);
    assert(ELEM_DIM < 4 && ELEM_DIM > 0);
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
    assert(ELEM_DIM < 4 && ELEM_DIM > 0);
    for(int i=0;i<ELEM_DIM+1;i++)
    {
        basisGradValues.push_back(ComputeBasisFunctionDerivative(psi,i));
    }

    return basisGradValues;    
}
