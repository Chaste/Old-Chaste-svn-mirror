#ifndef _TESTLINEARBASISFUNCTION_HPP_
#define _TESTLINEARBASISFUNCTION_HPP_

#include <cxxtest/TestSuite.h>
#include "LinearBasisFunction.cpp"
#include "GaussianQuadratureRule.hpp"
#include "VectorDouble.hpp"
#include <vector>

class TestLinearBasisFunction : public CxxTest::TestSuite 
{
	public:
	
	void testLinearBasisFunction0d() 
	{
		Point<0> zero;
		LinearBasisFunction<0> basis_func;
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 0), 1.0, 1e-12);

		std::vector<double> basis_function_vector;
		basis_function_vector = basis_func.ComputeBasisFunctions(zero);
		TS_ASSERT_EQUALS(basis_function_vector.size(), 1);
		TS_ASSERT_DELTA(basis_function_vector[0], 1.0, 1e-12);
		
		// check link with 0d quad rule works ok 
		GaussianQuadratureRule<0>  quad_rule(1);
		Point<0>   quad_point = quad_rule.GetQuadPoint(0);

		std::vector<double> basis_function_vector2;
		basis_function_vector2 = basis_func.ComputeBasisFunctions(quad_point);
		TS_ASSERT_EQUALS(basis_function_vector.size(), 1);
		TS_ASSERT_DELTA(basis_function_vector[0], 1.0, 1e-12);	
	}
	
	void testLinearBasisFunction1d()
	{
		Point<1> zero(0);
		Point<1> one(1);
		LinearBasisFunction<1> basis_func;
		
		// Single compute
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 0), 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 1), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(one,  0), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(one,  1), 1.0, 1e-12);
		
		// Mass compute
		std::vector<double> basis_function_vector;
		basis_function_vector = basis_func.ComputeBasisFunctions(zero);
		TS_ASSERT_EQUALS(basis_function_vector.size(), 2);
		TS_ASSERT_DELTA(basis_function_vector[0], 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_function_vector[1], 0.0, 1e-12);
		
		// Derivatives
		std::vector<VectorDouble> derivatives;
		derivatives = basis_func.ComputeBasisFunctionDerivatives(one);
		TS_ASSERT_EQUALS(derivatives.size(), 2);
		TS_ASSERT_DELTA(derivatives[0](0), -1, 1e-12);
		TS_ASSERT_DELTA(derivatives[1](0),  1, 1e-12);
	}
	
	void testLinearBasisFunction2d()
	{
		Point<2> zero(0,0);
		Point<2> onezero(1,0);
		Point<2> zeroone(0,1);
		
		LinearBasisFunction<2> basis_func;
		
		// Single compute
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 0), 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 1), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 2), 0.0, 1e-12);
		
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(onezero, 0), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(onezero, 1), 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(onezero, 2), 0.0, 1e-12);
		
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zeroone, 0), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zeroone, 1), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zeroone, 2), 1.0, 1e-12);
		
		// Mass compute
		std::vector<double> basis_function_vector;
		basis_function_vector = basis_func.ComputeBasisFunctions(zero);
		TS_ASSERT_EQUALS(basis_function_vector.size(), 3);
		TS_ASSERT_DELTA(basis_function_vector[0], 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_function_vector[1], 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_function_vector[2], 0.0, 1e-12);
		// Derivatives
		std::vector<VectorDouble> derivatives;
		derivatives = basis_func.ComputeBasisFunctionDerivatives(onezero);
		TS_ASSERT_EQUALS(derivatives.size(), 3);
		TS_ASSERT_DELTA(derivatives[0](0), -1, 1e-12);
		TS_ASSERT_DELTA(derivatives[1](0),  1, 1e-12);
		TS_ASSERT_DELTA(derivatives[2](0),  0, 1e-12);
	}


	void testLinearBasisFunction3d()
	{
		Point<3> zero(0,0,0);
		Point<3> zerozeroone(0,0,1);
		Point<3> zeroonezero(0,1,0);
		Point<3> onezerozero(1,0,0);
		
		LinearBasisFunction<3> basis_func;
		
		// Single compute
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 0), 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 1), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 2), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 3), 0.0, 1e-12);
		
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(onezerozero, 0), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(onezerozero, 1), 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(onezerozero, 2), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(onezerozero, 3), 0.0, 1e-12);
		
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zeroonezero, 0), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zeroonezero, 1), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zeroonezero, 2), 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(onezerozero, 3), 0.0, 1e-12);

		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerozeroone, 0), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerozeroone, 1), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerozeroone, 2), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerozeroone, 3), 1.0, 1e-12);
		
		// Mass compute
		std::vector<double> basis_function_vector;
		basis_function_vector = basis_func.ComputeBasisFunctions(zero);
		TS_ASSERT_EQUALS(basis_function_vector.size(), 4);
		TS_ASSERT_DELTA(basis_function_vector[0], 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_function_vector[1], 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_function_vector[2], 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_function_vector[3], 0.0, 1e-12);
		
		// Derivatives
		std::vector<VectorDouble> derivatives;
		derivatives = basis_func.ComputeBasisFunctionDerivatives(onezerozero);
		TS_ASSERT_EQUALS(derivatives.size(), 4);
		TS_ASSERT_DELTA(derivatives[0](0), -1, 1e-12);
		TS_ASSERT_DELTA(derivatives[1](0),  1, 1e-12);
		TS_ASSERT_DELTA(derivatives[2](0),  0, 1e-12);
		TS_ASSERT_DELTA(derivatives[3](0),  0, 1e-12);
	}
	
	
	void testComputeTransformedBasisFunctionDerivatives( void )
	{
		// 1D
		LinearBasisFunction<1> basis_func;
		Point<1> one(1);
		
		MatrixDouble invJ(1,1);
		invJ(0,0)=0.5;

		std::vector<VectorDouble> transDeriv = basis_func.ComputeTransformedBasisFunctionDerivatives(one,invJ);

		TS_ASSERT_EQUALS(transDeriv.size(), 2);
		TS_ASSERT_DELTA(transDeriv[0](0), -0.5, 1e-12);
		TS_ASSERT_DELTA(transDeriv[1](0),  0.5, 1e-12);
		
		// 2D
		LinearBasisFunction<2> basis_func2;
		Point<2> oneone(1,1);
		
		MatrixDouble invJ2(2,2);
		invJ2(0,0)=0.5;
		invJ2(1,1)=0.5;

		std::vector<VectorDouble> transDeriv2 = basis_func2.ComputeTransformedBasisFunctionDerivatives(oneone,invJ2);
		
		TS_ASSERT_EQUALS(transDeriv2.size(), 3);
		TS_ASSERT_DELTA(transDeriv2[0](0), -0.5, 1e-12);
		TS_ASSERT_DELTA(transDeriv2[1](0),  0.5, 1e-12);
		TS_ASSERT_DELTA(transDeriv2[2](0),    0, 1e-12);

		
		//3D
		LinearBasisFunction<3> basis_func3;
		Point<3> oneoneone(1,1,1);
		
		MatrixDouble invJ3(3,3);
		invJ3(0,0)=0.5;
		invJ3(1,1)=0.5;
		invJ3(2,2)=0.5;
		
		std::vector<VectorDouble> transDeriv3 = basis_func3.ComputeTransformedBasisFunctionDerivatives(oneoneone,invJ3);

		TS_ASSERT_EQUALS(transDeriv3.size(), 4);
		TS_ASSERT_DELTA(transDeriv3[0](0), -0.5, 1e-12);
		TS_ASSERT_DELTA(transDeriv3[1](0),  0.5, 1e-12);
		TS_ASSERT_DELTA(transDeriv3[2](0),    0, 1e-12);
		TS_ASSERT_DELTA(transDeriv3[3](0),    0, 1e-12);
			//TS_TRACE("here lin basis\n");	
	}
	
	void testComputeTransformedBasisFunction2( void )		
	{		
		// 2D - with better test data
		
		std::vector<Node<2>*> nodes;
		nodes.push_back(new Node<2>(0, false, 4.0, 3.0));
		nodes.push_back(new Node<2>(1, false, 6.0, 4.0));
		nodes.push_back(new Node<2>(2, false, 3.0, 5.0));
		Element<2,2> element(nodes);
		LinearBasisFunction<2> basis_function;
		
		const MatrixDouble *inverseJacobian = element.GetInverseJacobian();
		Point<2> evaluation_point(1,1); 
		std::vector<VectorDouble> trans_deriv =
			basis_function.ComputeTransformedBasisFunctionDerivatives(evaluation_point,
																		*inverseJacobian);
		
		TS_ASSERT_DELTA(trans_deriv[0](0),-0.2, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[0](1),-0.6, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[1](0),0.4, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[1](1),0.2, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[2](0),-0.2, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[2](1),0.4, 1e-12);
		
	}
};

#endif //_TESTLINEARBASISFUNCTION_HPP_
