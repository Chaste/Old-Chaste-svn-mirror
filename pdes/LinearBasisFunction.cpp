#ifndef _LINEARBASISFUNCTION_CPP_
#define _LINEARBASISFUNCTION_CPP_


#include "LinearBasisFunction.hpp"
#include "Point.hpp"
#include "MatrixDouble.hpp"
#include <cassert>


/**
 * Compute a basis function at a point within an element.
 * 
 * @param point The point at which to compute the basis function. The results
 *     are undefined if this is not within the canonical element.
 * @param basisIndex Which basis function to compute. This is a local index
 *     within a canonical element.
 * @return The value of the basis function.
 */
template <int ELEM_DIM>
double LinearBasisFunction<ELEM_DIM>::ComputeBasisFunction(Point<ELEM_DIM> point, int basisIndex) const
{
	assert(ELEM_DIM < 4 && ELEM_DIM >= 0);
	switch(ELEM_DIM)
	{
	case 0:
		assert(basisIndex == 0);
		return 1.0;
		break;			

	case 1:
		switch (basisIndex)
		{
			case 0:
				return 1.0 - point[0];
				break;
    		case 1:
        		return point[0];
        		break;
        	default:
        		assert(false);   
    	}
    	break;
    	
    case 2:
    	switch (basisIndex)
    	{
    		case 0:
    			return 1.0 - point[0] - point[1];
    			break;
    		case 1:
    			return point[0];
    			break;
    		case 2:
    			return point[1];
    			break;
    		default:
    			assert(false);   
    	}
    	break;
    	
    case 3:
    	switch (basisIndex)
    	{
    		case 0:	
        		return 1.0 - point[0] - point[1] - point[2];
        		break;
    		case 1:
    			return point[0];
    			break;
			case 2:
				return point[1];
				break;
			case 3:
				return point[2]; 
				break;
			default:
				assert(false);  
    	}
    	break;
	}
}


/**
 * Compute the derivative of a basis function at a point within an canonical element.
 * 
 * @param point The point at which to compute the basis function. The results
 *     are undefined if this is not within the canonical element.
 * @param basisIndex Which basis function to compute. This is a local index
 *     within a canonical element.
 * @return The derivative of the basis function. This is a vector (VectorDouble
 *     instance) giving the derivative along each axis.
 */
template <int ELEM_DIM>
VectorDouble LinearBasisFunction<ELEM_DIM>::ComputeBasisFunctionDerivative(Point<ELEM_DIM> point, int basisIndex) const
{
    VectorDouble gradN(ELEM_DIM);
    assert(ELEM_DIM < 4 && ELEM_DIM > 0);
    
	switch(ELEM_DIM)
	{
	case 1:
    	switch (basisIndex)
    	{
    		case 0:
    			gradN(0) =  - 1;
    			break;
			case 1:
				gradN(0) = 1; 
				break;
			default:
				assert(false);
    	}
    	break;
    	
    case 2:
    	switch (basisIndex)
    	{
    		case 0:
        		gradN(0) =  - 1;
        		gradN(1) = - 1;
        		break;
    		case 1:
        		gradN(0) = 1;
        		gradN(1) = 0;
        		break;
    		case 2:
    			gradN(0) = 0;
        		gradN(1) = 1;   
        		break;
    		default:
    			assert(false);
    	}
    	break;
    	
    case 3:
    	switch (basisIndex)
    	{
    		case 0:
        		gradN(0) =  - 1;
        		gradN(1) = - 1;
        		gradN(2) = - 1;
        		break;
    		case 1:
    			gradN(0) =  1;
	        	gradN(1) =  0;
    	    	gradN(2) =  0;
	    		break;
    		case 2:
    			gradN(0) =  0;
	        	gradN(1) =  1;
    	    	gradN(2) =  0;
    			break;
			case 3:
				gradN(0) =  0;
        		gradN(1) =  0;
        		gradN(2) =  1;
        		break;
    		default:
    			assert(false);   
    	}
    	break;
	}    
    return gradN;
}


/**
 * Compute all basis functions at a point within an element.
 * 
 * @param point The point at which to compute the basis functions. The results
 *     are undefined if this is not within the canonical element.
 * @return The values of the basis functions, in local index order.
 */
template <int ELEM_DIM>
std::vector<double> LinearBasisFunction<ELEM_DIM>::ComputeBasisFunctions(Point<ELEM_DIM> point) const
{
    std::vector<double> basisValues(ELEM_DIM+1);
    assert(ELEM_DIM < 4 && ELEM_DIM >= 0);
    for(int i=0;i<ELEM_DIM+1;i++)
    {
        basisValues[i] = ComputeBasisFunction(point,i);
    }
    return basisValues;
}



/**
 * Compute the derivatives of all basis functions at a point within an element.
 * 
 * @param point The point at which to compute the basis functions. The results
 *     are undefined if this is not within the canonical element.
 * @return The derivatives of the basis functions, in local index order. Each
 *     entry is a vector (VectorDouble instance) giving the derivative along
 *     each axis.
 */
template <int ELEM_DIM>
std::vector<VectorDouble>  LinearBasisFunction<ELEM_DIM>::ComputeBasisFunctionDerivatives(Point<ELEM_DIM> point) const
{
    std::vector<VectorDouble> basisGradValues;
    assert(ELEM_DIM < 4 && ELEM_DIM > 0);
    for(int i=0;i<ELEM_DIM+1;i++)
    {
        basisGradValues.push_back(ComputeBasisFunctionDerivative(point,i));
    }

    return basisGradValues;    
}


/**
 * Compute the derivatives of all basis functions at a point within an element.
 * This method will transform the results, for use within gaussian quadrature
 * for example.
 * 
 * @param point The point at which to compute the basis functions. The results
 *     are undefined if this is not within the canonical element.
 * @param inverseJacobian The inverse of the Jacobian matrix mapping the real
 *     element into the canonical element.
 * @return The derivatives of the basis functions, in local index order. Each
 *     entry is a vector (VectorDouble instance) giving the derivative along
 *     each axis.
 */
template <int ELEM_DIM>
std::vector<VectorDouble>  LinearBasisFunction<ELEM_DIM>::ComputeTransformedBasisFunctionDerivatives(Point<ELEM_DIM> point, MatrixDouble inverseJacobian) const
{
    assert(ELEM_DIM < 4 && ELEM_DIM > 0);
    assert(inverseJacobian.Rows()==inverseJacobian.Columns());
    std::vector<VectorDouble> basisGradValues = ComputeBasisFunctionDerivatives(point);
  	std::vector<VectorDouble> transformedGradValues;
    
    for(int i=0;i<ELEM_DIM+1;i++)
    {
        transformedGradValues.push_back( basisGradValues[i]*inverseJacobian );
    }

    return transformedGradValues;    	
}


#endif // _LINEARBASISFUNCTION_CPP_
