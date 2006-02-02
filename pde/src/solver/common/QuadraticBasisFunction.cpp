#ifndef _QUADRATICBASISFUNCTION_CPP_
#define _QUADRATICBASISFUNCTION_CPP_

#include "QuadraticBasisFunction.hpp"
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
double QuadraticBasisFunction<ELEM_DIM>::ComputeBasisFunction(const Point<ELEM_DIM> &rPoint, int basisIndex) const
{
	assert(ELEM_DIM < 4 && ELEM_DIM >= 0);
	double x, y, z;
	switch(ELEM_DIM)
	{
	case 0:
		assert(basisIndex == 0);
		return 1.0;
		break;			

	case 1:
		x = rPoint[0];
		switch (basisIndex)
		{
			case 0:
        		return 2.0*(x-1.0)*(x-0.5);
        		break;
    		case 1:
        		return 2.0*x*(x-0.5);
        		break;
    		case 2:
        		return 4.0*x*(1.0-x);
        		break;
    		default:
    			assert(false);
		}
    	break;
    	
    case 2:
    	x = rPoint[0];
    	y = rPoint[1];
    	switch (basisIndex)
    	{
    		case 0:
         		return 2.0 * (1.0 - x - y) * (0.5 - x - y);
         		break;
    	  	case 1:
        		return 2.0*x*(x-0.5);
        		break;
    		case 2:
        		return 2.0*y*(y-0.5);
        		break;
    		case 3:
        		return 4.0 * (1.0 - x - y) * x;
        		break;
    		case 4:
        		return 4.0 * (1.0 - x - y) * y;
        		break;
    		case 5:
        		return 4.0 * y * x;
        		break;
    		default:
    			assert(false);
    	}
    	break;
    	
    case 3:
    	x = rPoint[0];
    	y = rPoint[1];
    	z = rPoint[2];
        switch (basisIndex)
        {
        	case 0:
    			return 2.0 * (1.0 - x - y - z) * (0.5 - x - y - z);
    			break;
    		case 1:
        		return 2.0*x*(x-0.5);
        		break;
    		case 2:
        		return 2.0*y*(y-0.5);
        		break;
    		case 3:
        		return 2.0*z*(z-0.5);
        		break;
    		case 4:
    			return 4.0 * (1.0 - x - y - z) * x;
    			break;
    		case 5:
        		return 4.0 * (1.0 - x - y - z) * y;
        		break;
    		case 6:
        		return 4.0 * (1.0 - x - y - z) * z;
        		break;
        	case 7:
        		return 4.0 * y * x;
        		break;
    		case 8:
        		return 4.0 * x * z;
        		break;
    		case 9:
    			return 4.0 * y * z;
    			break;
    		default:
    			assert(false);
        }
    	break;
	}
	return 0.0; // Avoid compiler warning
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
VectorDouble QuadraticBasisFunction<ELEM_DIM>::ComputeBasisFunctionDerivative(const Point<ELEM_DIM> &rPoint, int basisIndex) const
{
    VectorDouble gradN(ELEM_DIM);
    assert(ELEM_DIM < 4 && ELEM_DIM > 0);
    
    double x, y, z;
	switch(ELEM_DIM)
	{
	case 1:
		x=rPoint[0];
		switch (basisIndex)
		{
			case 0:
        		gradN(0) =  4.0*x-3.0;
        		break;
    		case 1:
    			gradN(0) = 4.0*x-1.0;
    			break;
    		case 2:
    			gradN(0) = 4.0-8.0*x;
    			break;
    		default:
    			assert(false);
		}
    	break;
    	
    case 2:
    	x = rPoint[0];
    	y = rPoint[1];
        switch (basisIndex)
        {
        	case 0:
    			gradN(0) = -3.0 + 4.0*x + 4.0*y;
        		gradN(1) = -3.0 + 4.0*x + 4.0*y;
        		break;
        	case 1:
        		gradN(0) = 4.0*x - 1.0;
        		gradN(1) = 0.0;
        		break;
        	case 2:
        		gradN(0) = 0.0;
        		gradN(1) = 4.0*y - 1.0;
        		break;
        	case 3:
        		gradN(0) = 4.0-8.0*x-4.0*y;
        		gradN(1) = -4.0*x;
        		break;
        	case 4:
    			gradN(0) = -4.0*y;
        		gradN(1) = 4.0-4.0*x-8.0*y;
        		break;
        	case 5:
    			gradN(0) = 4.0*y;
        		gradN(1) = 4.0*x;
        		break;
        	default:
        		assert (false);
    	}
    	break;
    	
    case 3:
    	x = rPoint[0];
    	y = rPoint[1];
    	z = rPoint[2];
    	switch (basisIndex)
    	{
    		case 0:
    		   	gradN(0) = -3.0 + 4.0*(x+y+z);
        		gradN(1) = -3.0 + 4.0*(x+y+z);
        		gradN(2) = -3.0 + 4.0*(x+y+z);
        		break;
        	case 1:
        		gradN(0) =  4.0*x-1.0;
        		gradN(1) =  0;
        		gradN(2) =  0;
        		break;
        	case 2:
        		gradN(0) =  0;
        		gradN(1) =  4.0*y-1.0;
        		gradN(2) =  0;
        		break;
        	case 3:
        		gradN(0) =  0;
        		gradN(1) =  0;
        		gradN(2) =  4.0*z-1.0;
        		break;
        	case 4:
        		gradN(0) =  4.0-8.0*x-4.0*y-4.0*z;
        		gradN(1) =  -4.0*x;
        		gradN(2) =  -4.0*x;
        		break;
        	case 5:
        		gradN(0) =  -4.0*y;
        		gradN(1) =  4.0-4.0*x-8.0*y-4.0*z;
        		gradN(2) =  -4.0*y; 
        		break;
        	case 6: 
        		gradN(0) =  -4.0*z;
        		gradN(1) =  -4.0*z;
        		gradN(2) =  4.0-4.0*x-4.0*y-8.0*z;
        		break;
        	case 7:
        		gradN(0) =  4.0*y;
        		gradN(1) =  4.0*x;
        		gradN(2) =  0.0;
        		break;
        	case 8:   
    			gradN(0) =  4.0*z;
        		gradN(1) =  0;
        		gradN(2) =  4.0*x;
        		break;
        	case 9:
        		gradN(0) =  0;
        		gradN(1) =  4.0*z;
        		gradN(2) =  4.0*y;
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
std::vector<double> QuadraticBasisFunction<ELEM_DIM>::ComputeBasisFunctions(const Point<ELEM_DIM> &rPoint) const
{
    assert(ELEM_DIM < 4 && ELEM_DIM >= 0);
    std::vector<double> basisValues((ELEM_DIM+1)*(ELEM_DIM+2)/2);
    for(int i=0;i<(ELEM_DIM+1)*(ELEM_DIM+2)/2;i++)
    {
        basisValues[i] = ComputeBasisFunction(rPoint,i);
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
std::vector<VectorDouble>  QuadraticBasisFunction<ELEM_DIM>::ComputeBasisFunctionDerivatives(const Point<ELEM_DIM> &rPoint) const
{
    assert(ELEM_DIM < 4 && ELEM_DIM > 0);
    std::vector<VectorDouble> basisGradValues;
    basisGradValues.reserve((ELEM_DIM+1)*(ELEM_DIM+2)/2);
    for(int i=0; i<(ELEM_DIM+1)*(ELEM_DIM+2)/2; i++)
    {
        basisGradValues.push_back(ComputeBasisFunctionDerivative(rPoint, i));
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
std::vector<VectorDouble> QuadraticBasisFunction<ELEM_DIM>::ComputeTransformedBasisFunctionDerivatives(const Point<ELEM_DIM> &rPoint, const MatrixDouble &rInverseJacobian) const
{
    assert(ELEM_DIM < 4 && ELEM_DIM > 0);
    assert(rInverseJacobian.Rows() == rInverseJacobian.Columns());
    std::vector<VectorDouble> basisGradValues = ComputeBasisFunctionDerivatives(rPoint);
  	std::vector<VectorDouble> transformedGradValues;
  	transformedGradValues.reserve((ELEM_DIM+1)*(ELEM_DIM+2)/2);
    
    for(int i=0; i<(ELEM_DIM+1)*(ELEM_DIM+2)/2; i++)
    {
        transformedGradValues.push_back(basisGradValues[i] * rInverseJacobian);
    }

    return transformedGradValues;    	
}

#endif // _QUADRATICBASISFUNCTION_CPP_


