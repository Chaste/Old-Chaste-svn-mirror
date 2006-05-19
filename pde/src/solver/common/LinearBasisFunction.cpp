#ifndef _LINEARBASISFUNCTION_CPP_
#define _LINEARBASISFUNCTION_CPP_


#include "LinearBasisFunction.hpp"
#include "Point.hpp"
#include "MatrixDouble.hpp"
#include "VectorDoubleUblasConverter.hpp"
#include "MatrixDoubleUblasConverter.hpp"
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
double LinearBasisFunction<ELEM_DIM>::ComputeBasisFunction(const Point<ELEM_DIM> &rPoint, int basisIndex) const
{
	assert(ELEM_DIM < 4 && ELEM_DIM > 0);
	switch(ELEM_DIM)
	{
	case 1:
		switch (basisIndex)
		{
			case 0:
				return 1.0 - rPoint[0];
				break;
    		case 1:
        		return rPoint[0];
        		break;
        	default:
        		assert(false);   
    	}
    	break;
    	
    case 2:
    	switch (basisIndex)
    	{
    		case 0:
    			return 1.0 - rPoint[0] - rPoint[1];
    			break;
    		case 1:
    			return rPoint[0];
    			break;
    		case 2:
    			return rPoint[1];
    			break;
    		default:
    			assert(false);   
    	}
    	break;
    	
    case 3:
    	switch (basisIndex)
    	{
    		case 0:	
        		return 1.0 - rPoint[0] - rPoint[1] - rPoint[2];
        		break;
    		case 1:
    			return rPoint[0];
    			break;
			case 2:
				return rPoint[1];
				break;
			case 3:
				return rPoint[2]; 
				break;
			default:
				assert(false);  
    	}
    	break;
	}
	return 0.0; // Avoid compiler warning
}

/**
 * Compute a basis function at a point within an element.
 * 
 * @param point The point at which to compute the basis function. The results
 *     are undefined if this is not within the canonical element.
 * @param basisIndex Which basis function to compute. This is a local index
 *     within a canonical element.
 * @return The value of the basis function.
 */
double LinearBasisFunction<0>::ComputeBasisFunction(const Point<0> &rPoint, int basisIndex) const
{
    assert(basisIndex == 0);
    return 1.0;
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
c_vector<double, ELEM_DIM> LinearBasisFunction<ELEM_DIM>::ComputeBasisFunctionDerivative(const Point<ELEM_DIM>&, int basisIndex) const
{
    c_vector<double, ELEM_DIM> gradN;
    assert(ELEM_DIM < 4 && ELEM_DIM > 0);
    
	switch(ELEM_DIM)
	{
	case 1:
    	switch (basisIndex)
    	{
    		case 0:
    			gradN(0) = -1;
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
        		gradN(0) = -1;
        		gradN(1) = -1;
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
        		gradN(0) = -1;
        		gradN(1) = -1;
        		gradN(2) = -1;
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
c_vector<double, ELEM_DIM+1> LinearBasisFunction<ELEM_DIM>::ComputeBasisFunctions(const Point<ELEM_DIM> &rPoint) const
{
    c_vector<double, ELEM_DIM+1> basisValues;
    assert(ELEM_DIM < 4 && ELEM_DIM >= 0);
    for(int i=0;i<ELEM_DIM+1;i++)
    {
        basisValues(i) = ComputeBasisFunction(rPoint, i);
    }
    return basisValues;
}

/**
 * Compute all basis functions at a point within an element.
 * 
 * @param point The point at which to compute the basis functions. The results
 *     are undefined if this is not within the canonical element.
 * @return The values of the basis functions, in local index order.
 */
c_vector<double, 1> LinearBasisFunction<0>::ComputeBasisFunctions(const Point<0> &rPoint) const
{
    c_vector<double, 1> basisValues;
    basisValues(0) = ComputeBasisFunction(rPoint, 0);
    return basisValues;
}

/**
 * Compute the derivatives of all basis functions at a point within an element.
 * 
 * @param point The point at which to compute the basis functions. The results
 *     are undefined if this is not within the canonical element.
 * @return The derivatives of the basis functions as the column vectors of a matrix
 *     in local index order. 
 */
template <int ELEM_DIM>
c_matrix<double, ELEM_DIM, ELEM_DIM+1>  LinearBasisFunction<ELEM_DIM>::ComputeBasisFunctionDerivatives(const Point<ELEM_DIM> &rPoint) const
{
    assert(ELEM_DIM < 4 && ELEM_DIM > 0);
    c_matrix<double, ELEM_DIM, ELEM_DIM+1> basisGradValues;
    
    
    for(unsigned j=0;j<ELEM_DIM+1;j++)
    {
        matrix_column< c_matrix<double, ELEM_DIM, ELEM_DIM+1> >  column(basisGradValues,j);
        column=ComputeBasisFunctionDerivative(rPoint,j);
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
c_matrix<double, ELEM_DIM, ELEM_DIM+1> LinearBasisFunction<ELEM_DIM>::ComputeTransformedBasisFunctionDerivatives(const Point<ELEM_DIM> &rPoint, const c_matrix<double, ELEM_DIM, ELEM_DIM> &rInverseJacobian) const
{
    assert(ELEM_DIM < 4 && ELEM_DIM > 0);
    c_matrix<double, ELEM_DIM, ELEM_DIM+1> basisGradValues = ComputeBasisFunctionDerivatives(rPoint);
  	std::vector<c_vector<double, ELEM_DIM> > transformedGradValues;
    
    return prod(trans(rInverseJacobian), basisGradValues);    	
}


#endif // _LINEARBASISFUNCTION_CPP_
