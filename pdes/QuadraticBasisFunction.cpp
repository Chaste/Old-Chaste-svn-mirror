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
double QuadraticBasisFunction<ELEM_DIM>::ComputeBasisFunction(Point<ELEM_DIM> point, int basisIndex)
{
	assert(ELEM_DIM < 4 && ELEM_DIM >= 0);
	switch(ELEM_DIM)
	{
	case 0:
		assert(basisIndex == 0);
		return 1.0;
		break;			

	case 1:
		assert(basisIndex == 0 || basisIndex == 1 || basisIndex == 2);
    	if(basisIndex == 0)
    	{	
    		double x = point[0];
        	return 2.0*(x-1.0)*(x-0.5);
    	}
    	else if(basisIndex == 1)
    	{
    		double x = point[0];
        	return 2.0*x*(x-0.5);
    	}
    	else
    	{
        	double x = point[0];
        	return 4.0*x*(1.0-x);
    	}
    	break;
    	
    case 2:
    	assert(basisIndex == 0 || basisIndex == 1 || basisIndex == 2 || basisIndex == 3 || basisIndex == 4 || basisIndex == 5);
    	if(basisIndex == 0)
    	{	
    		double x = point[0];
    		double y = point[1];
        	return 2.0 * (1.0 - x - y) * (0.5 - x - y);
    	}
    	else if(basisIndex == 1)
    	{	
        	double x = point[0];
        	return 2.0*x*(x-0.5);
    	}
    	else if(basisIndex == 2)
    	{
    		double y = point[1];
        	return 2.0*y*(y-0.5);
    	}
    	else if(basisIndex == 3)
    	{	
    		double x = point[0];
    		double y = point[1];
        	return 4.0 * (1.0 - x - y) * x;
    	}
    	else if(basisIndex == 4)
    	{	
    		double x = point[0];
    		double y = point[1];
        	return 4.0 * y * x;
    	}
    	else
    	{
    		double x = point[0];
    		double y = point[1];
        	return 4.0 * (1.0 - x - y) * y;
    	}
    	break;
    	
    case 3:
    	assert(basisIndex == 0 || basisIndex == 1 || basisIndex == 2 || basisIndex == 3 || basisIndex == 4 || basisIndex == 5 || basisIndex == 6 || basisIndex == 7 || basisIndex == 8 || basisIndex == 9);
    	if(basisIndex == 0)
    	{	
    		double x = point[0];
    		double y = point[1];
    		double z = point[2];
        	return 2.0 * (1.0 - x - y - z) * (0.5 - x - y - z);
    	}
    	else if(basisIndex == 1)
    	{	
    		double x = point[0];
        	return 2.0*x*(x-0.5);
    	}
    	else if(basisIndex == 2)
    	{	
        	double y = point[1];
        	return 2.0*y*(y-0.5);
    	}
    	else if(basisIndex == 3)
    	{
        	double z = point[2];
        	return 2.0*z*(z-0.5);
    	}
    	else if(basisIndex == 4)
    	{	
    		double x = point[0];
    		double y = point[1];
    		double z = point[2];
        	return 4.0 * (1.0 - x - y - z) * x;
    	}
    	else if(basisIndex == 5)
    	{	
    		double x = point[0];
    		double y = point[1];
        	return 4.0 * y * x;
    	}
    	else if(basisIndex == 6)
    	{	
    		double x = point[0];
    		double y = point[1];
    		double z = point[2];
        	return 4.0 * (1.0 - x - y - z) * y;
    	}
    	else if(basisIndex == 7)
    	{	
    		double x = point[0];
    		double y = point[1];
    		double z = point[2];
        	return 4.0 * (1.0 - x - y - z) * z;
    	}
    	else if(basisIndex == 8)
    	{	
    		double x = point[0];
    		double z = point[2];
        	return 4.0 * x * z;
    	}
    	else
    	{	
    		double y = point[1];
    		double z = point[2];
        	return 4.0 * y * z;
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
VectorDouble QuadraticBasisFunction<ELEM_DIM>::ComputeBasisFunctionDerivative(Point<ELEM_DIM> point, int basisIndex)
{
    VectorDouble gradN(ELEM_DIM);
    assert(ELEM_DIM < 4 && ELEM_DIM > 0);
    
	switch(ELEM_DIM)
	{
	case 1:
		assert(basisIndex == 0 || basisIndex == 1 || basisIndex == 2);
    	if(basisIndex == 0)
    	{	double x = point[0];
        	gradN(0) =  4.0*x-3.0;
    	}
    	else if(basisIndex == 1)
    	{
    		double x = point[0];
        	gradN(0) = 4.0*x-1.0;   
    	}
    	else
    	{
    		double x = point[0];
    		gradN(0) = 4.0-8.0*x;
    	}
    	break;
    	
    case 2:
    	assert(basisIndex == 0 || basisIndex == 1 || basisIndex == 2 || basisIndex == 3 || basisIndex == 4 || basisIndex == 5);
    	if(basisIndex == 0)
    	{	
    		double x = point[0];
    		double y = point[1];
        	gradN(0) = -3.0 + 4.0*x + 4.0*y;
        	gradN(1) = -3.0 + 4.0*x + 4.0*y;
    	}
    	else if(basisIndex == 1)
    	{	
    		double x = point[0];
        	gradN(0) = 4.0*x - 1.0;
        	gradN(1) = 0.0;
    	}
    	else if(basisIndex == 2)
    	{
    		double y = point[1];
        	gradN(0) = 0.0;
        	gradN(1) = 4.0*y - 1.0;
    	}
    	else if(basisIndex == 3)
    	{
   			double x = point[0];
    		double y = point[1];
        	gradN(0) = 4.0-8.0*x-4.0*y;
        	gradN(1) = -4.0*x;
    	}
    	else if(basisIndex == 4)
    	{
   			double x = point[0];
    		double y = point[1];
        	gradN(0) = 4.0*y;
        	gradN(1) = 4.0*x;
    	}
    	else
    	{
   			double x = point[0];
    		double y = point[1];
        	gradN(0) = -4.0*y;
        	gradN(1) = 4.0-4.0*x-8.0*y;
    	}
    	break;
    	
    case 3:
    	assert(basisIndex == 0 || basisIndex == 1 || basisIndex == 2 || basisIndex == 3 || basisIndex == 4 || basisIndex == 5 || basisIndex == 6 || basisIndex == 7 || basisIndex == 8 || basisIndex == 9);
    	if(basisIndex == 0)
    	{	
    		double x = point[0];
    		double y = point[1];
    		double z = point[2];
        	gradN(0) = -3.0 + 4.0*(x+y+z);
        	gradN(1) = -3.0 + 4.0*(x+y+z);
        	gradN(2) = -3.0 + 4.0*(x+y+z);
    	}
    	else if(basisIndex == 1)
    	{	
    		double x = point[0];
    		double y = point[1];
    		double z = point[2];
        	gradN(0) =  4.0*x-1.0;
        	gradN(1) =  0;
        	gradN(2) =  0;
    	}
    	else if(basisIndex == 2)
    	{	
    		double x = point[0];
    		double y = point[1];
    		double z = point[2];
        	gradN(0) =  0;
        	gradN(1) =  4.0*y-1.0;
        	gradN(2) =  0;
    	}
    	else if(basisIndex == 3)
    	{
    		double x = point[0];
    		double y = point[1];
    		double z = point[2];
        	gradN(0) =  0;
        	gradN(1) =  0;
        	gradN(2) =  4.0*z-1.0;   
    	}
    	else if(basisIndex == 4)
    	{
    		double x = point[0];
    		double y = point[1];
    		double z = point[2];
        	gradN(0) =  4.0-8.0*x-4.0*y-4.0*z;
        	gradN(1) =  -4.0*x;
        	gradN(2) =  -4.0*x;   
    	}
    	else if(basisIndex == 5)
    	{
    		double x = point[0];
    		double y = point[1];
    		double z = point[2];
        	gradN(0) =  4.0*y;
        	gradN(1) =  4.0*x;
        	gradN(2) =  0.0;   
    	}
    	else if(basisIndex == 6)
    	{
    		double x = point[0];
    		double y = point[1];
    		double z = point[2];
        	gradN(0) =  -4.0*y;
        	gradN(1) =  4.0-4.0*x-8.0*y-4.0*z;
        	gradN(2) =  -4.0*y;   
    	}
    	else if(basisIndex == 7)
    	{
    		double x = point[0];
    		double y = point[1];
    		double z = point[2];
        	gradN(0) =  -4.0*z;
        	gradN(1) =  -4.0*z;
        	gradN(2) =  4.0-4.0*x-4.0*y-8.0*z;   
    	}
    	else if(basisIndex == 8)
    	{
    		double x = point[0];
    		double y = point[1];
    		double z = point[2];
        	gradN(0) =  4.0*z;
        	gradN(1) =  0;
        	gradN(2) =  4.0*x;   
    	}
    	else
    	{
    		double x = point[0];
    		double y = point[1];
    		double z = point[2];
        	gradN(0) =  0;
        	gradN(1) =  4.0*z;
        	gradN(2) =  4.0*y;   
    	}
    	
    	
    	
    	break;
	}    
    return gradN;
}
//
//
///**
// * Compute all basis functions at a point within an element.
// * 
// * @param point The point at which to compute the basis functions. The results
// *     are undefined if this is not within the canonical element.
// * @return The values of the basis functions, in local index order.
// */
template <int ELEM_DIM>
std::vector<double> QuadraticBasisFunction<ELEM_DIM>::ComputeBasisFunctions(Point<ELEM_DIM> point)
{
    std::vector<double> basisValues((ELEM_DIM+1)*(ELEM_DIM+2)/2);
    assert(ELEM_DIM < 4 && ELEM_DIM >= 0);
    for(int i=0;i<(ELEM_DIM+1)*(ELEM_DIM+2)/2;i++)
    {
        basisValues[i] = ComputeBasisFunction(point,i);
    }
    return basisValues;
}
//
//
//
///**
// * Compute the derivatives of all basis functions at a point within an element.
// * 
// * @param point The point at which to compute the basis functions. The results
// *     are undefined if this is not within the canonical element.
// * @return The derivatives of the basis functions, in local index order. Each
// *     entry is a vector (VectorDouble instance) giving the derivative along
// *     each axis.
// */
template <int ELEM_DIM>
std::vector<VectorDouble>  QuadraticBasisFunction<ELEM_DIM>::ComputeBasisFunctionDerivatives(Point<ELEM_DIM> point)
{
    std::vector<VectorDouble> basisGradValues;
    assert(ELEM_DIM < 4 && ELEM_DIM > 0);
    for(int i=0;i<(ELEM_DIM+1)*(ELEM_DIM+2)/2;i++)
    {
        basisGradValues.push_back(ComputeBasisFunctionDerivative(point,i));
    }

    return basisGradValues;    
}
//
//
///**
// * Compute the derivatives of all basis functions at a point within an element.
// * This method will transform the results, for use within gaussian quadrature
// * for example.
// * 
// * @param point The point at which to compute the basis functions. The results
// *     are undefined if this is not within the canonical element.
// * @param inverseJacobian The inverse of the Jacobian matrix mapping the real
// *     element into the canonical element.
// * @return The derivatives of the basis functions, in local index order. Each
// *     entry is a vector (VectorDouble instance) giving the derivative along
// *     each axis.
// */
template <int ELEM_DIM>
std::vector<VectorDouble>  QuadraticBasisFunction<ELEM_DIM>::ComputeTransformedBasisFunctionDerivatives(Point<ELEM_DIM> point, MatrixDouble inverseJacobian)
{
    assert(ELEM_DIM < 4 && ELEM_DIM > 0);
    assert(inverseJacobian.Rows()==inverseJacobian.Columns());
    std::vector<VectorDouble> basisGradValues = ComputeBasisFunctionDerivatives(point);
  	std::vector<VectorDouble> transformedGradValues;
    
    for(int i=0;i<(ELEM_DIM+1)*(ELEM_DIM+2)/2;i++)
    {
        transformedGradValues.push_back( basisGradValues[i]*inverseJacobian );
    }

    return transformedGradValues;    	
}

#endif // _QUADRATICBASISFUNCTION_CPP_


