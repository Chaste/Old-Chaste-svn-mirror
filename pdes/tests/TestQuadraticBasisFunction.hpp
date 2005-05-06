#ifndef _TESTQUADRATICBASISFUNCTION_HPP_
#define _TESTQUADRATICBASISFUNCTION_HPP_

#include <cxxtest/TestSuite.h>
#include "QuadraticBasisFunction.cpp"
#include "GaussianQuadratureRule.hpp"
#include "VectorDouble.hpp"
#include "BasisFunctionsCheckers.hpp"
#include <vector>


class TestQuadraticBasisFunction : public CxxTest::TestSuite 
{
	public:
	
	void TestQuadraticBasisFunction0d() 
	{
		Point<0> zero;
		QuadraticBasisFunction<0> basis_func;
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 0), 1.0, 1e-12);
	}
		
	void TestQuadraticBasisFunction1d()
	{
		std::vector<Point<1>*> evaluation_points;
		Point<1> zero(0);
		Point<1> one(1);
		Point<1> half(0.5);
		evaluation_points.push_back(&zero);
		evaluation_points.push_back(&one);
		evaluation_points.push_back(&half);
		
		QuadraticBasisFunction<1> basis_func;
		
		BasisFunctionsCheckers<1> checker;
		checker.checkBasisFunctions(&basis_func, evaluation_points);
		
		// Derivatives
		std::vector<VectorDouble> derivatives;
		derivatives = basis_func.ComputeBasisFunctionDerivatives(zero);
		TS_ASSERT_EQUALS(derivatives.size(), 3);
		TS_ASSERT_DELTA(derivatives[0](0), -3.0, 1e-12);
		TS_ASSERT_DELTA(derivatives[1](0),  -1.0, 1e-12);
		TS_ASSERT_DELTA(derivatives[2](0),  4.0, 1e-12);
		
		derivatives = basis_func.ComputeBasisFunctionDerivatives(one);
		TS_ASSERT_EQUALS(derivatives.size(), 3);
		TS_ASSERT_DELTA(derivatives[0](0), 1.0, 1e-12);
		TS_ASSERT_DELTA(derivatives[1](0),  3.0, 1e-12);
		TS_ASSERT_DELTA(derivatives[2](0),  -4.0, 1e-12);
		
		derivatives = basis_func.ComputeBasisFunctionDerivatives(half);
		TS_ASSERT_EQUALS(derivatives.size(), 3);
		TS_ASSERT_DELTA(derivatives[0](0), -1.0, 1e-12);
		TS_ASSERT_DELTA(derivatives[1](0),  1.0, 1e-12);
		TS_ASSERT_DELTA(derivatives[2](0),  0.0, 1e-12);
	}
	
	void TestQuadraticBasisFunction2d()
	{
		Point<2> zero(0,0);
		Point<2> onezero(1,0);
		Point<2> zeroone(0,1);
		Point<2> halfzero(0.5,0);
		Point<2> zerohalf(0,0.5);
		Point<2> halfhalf(0.5,0.5);
		
		std::vector<Point<2>*> evaluation_points;
		evaluation_points.push_back(&zero);
		evaluation_points.push_back(&onezero);
		evaluation_points.push_back(&zeroone);
		evaluation_points.push_back(&halfzero);
		evaluation_points.push_back(&zerohalf);
		evaluation_points.push_back(&halfhalf);
		
		QuadraticBasisFunction<2> basis_func;
		
		BasisFunctionsCheckers<2> checker;
		checker.checkBasisFunctions(&basis_func, evaluation_points);
				
		// Derivatives
		std::vector<VectorDouble> derivatives;
		derivatives = basis_func.ComputeBasisFunctionDerivatives(onezero);
		TS_ASSERT_EQUALS(derivatives.size(), 6);
		TS_ASSERT_DELTA(derivatives[0](0), 1, 1e-12);
		TS_ASSERT_DELTA(derivatives[1](0), 3, 1e-12);
		TS_ASSERT_DELTA(derivatives[2](0), 0, 1e-12);
		TS_ASSERT_DELTA(derivatives[3](0), -4, 1e-12);
		TS_ASSERT_DELTA(derivatives[4](0), 0, 1e-12);
		TS_ASSERT_DELTA(derivatives[5](0), 0, 1e-12);
		TS_ASSERT_DELTA(derivatives[0](1), 1, 1e-12);
		TS_ASSERT_DELTA(derivatives[1](1), 0, 1e-12);
		TS_ASSERT_DELTA(derivatives[2](1), -1, 1e-12);
		TS_ASSERT_DELTA(derivatives[3](1), -4, 1e-12);
		TS_ASSERT_DELTA(derivatives[4](1), 0, 1e-12);
		TS_ASSERT_DELTA(derivatives[5](1), 4, 1e-12);
		
		derivatives = basis_func.ComputeBasisFunctionDerivatives(zero);
		TS_ASSERT_EQUALS(derivatives.size(), 6);
		TS_ASSERT_DELTA(derivatives[0](0), -3, 1e-12);
		TS_ASSERT_DELTA(derivatives[1](0), -1, 1e-12);
		TS_ASSERT_DELTA(derivatives[2](0), 0, 1e-12);
		TS_ASSERT_DELTA(derivatives[3](0), 4, 1e-12);
		TS_ASSERT_DELTA(derivatives[4](0), 0, 1e-12);
		TS_ASSERT_DELTA(derivatives[5](0), 0, 1e-12);
		TS_ASSERT_DELTA(derivatives[0](1), -3, 1e-12);
		TS_ASSERT_DELTA(derivatives[1](1), 0, 1e-12);
		TS_ASSERT_DELTA(derivatives[2](1), -1, 1e-12);
		TS_ASSERT_DELTA(derivatives[3](1), 0, 1e-12);
		TS_ASSERT_DELTA(derivatives[4](1), 4, 1e-12);
		TS_ASSERT_DELTA(derivatives[5](1), 0, 1e-12);
		
		derivatives = basis_func.ComputeBasisFunctionDerivatives(zeroone);
		TS_ASSERT_EQUALS(derivatives.size(), 6);
		TS_ASSERT_DELTA(derivatives[0](0), 1, 1e-12);
		TS_ASSERT_DELTA(derivatives[1](0), -1, 1e-12);
		TS_ASSERT_DELTA(derivatives[2](0), 0, 1e-12);
		TS_ASSERT_DELTA(derivatives[3](0), 0, 1e-12);
		TS_ASSERT_DELTA(derivatives[4](0), -4, 1e-12);
		TS_ASSERT_DELTA(derivatives[5](0), 4, 1e-12);
		TS_ASSERT_DELTA(derivatives[0](1), 1, 1e-12);
		TS_ASSERT_DELTA(derivatives[1](1), 0, 1e-12);
		TS_ASSERT_DELTA(derivatives[2](1), 3, 1e-12);
		TS_ASSERT_DELTA(derivatives[3](1), 0, 1e-12);
		TS_ASSERT_DELTA(derivatives[4](1), -4, 1e-12);
		TS_ASSERT_DELTA(derivatives[5](1), 0, 1e-12);
	}

	void TestQuadraticBasisFunction3d()
	{
		Point<3> zero(0,0,0);
		Point<3> zerozeroone(0,0,1);
		Point<3> zeroonezero(0,1,0);
		Point<3> onezerozero(1,0,0);
		Point<3> halfhalfzero(0.5,0.5,0);
		Point<3> halfzerozero(0.5,0,0);
		Point<3> halfzerohalf(0.5,0.0,0.5);
		Point<3> zerohalfhalf(0.0,0.5,0.5);
		Point<3> zerohalfzero(0.0,0.5,0.0);
		Point<3> zerozerohalf(0.0,0.0,0.5);
		
		std::vector<Point<3>*> evaluation_points;
		evaluation_points.push_back(&zero);
		evaluation_points.push_back(&onezerozero);
		evaluation_points.push_back(&zeroonezero);
		evaluation_points.push_back(&zerozeroone);
		evaluation_points.push_back(&halfzerozero);
		evaluation_points.push_back(&zerohalfzero);
		evaluation_points.push_back(&zerozerohalf);
		evaluation_points.push_back(&halfhalfzero);
		evaluation_points.push_back(&halfzerohalf);
		evaluation_points.push_back(&zerohalfhalf);

		QuadraticBasisFunction<3> basis_func;
		
		BasisFunctionsCheckers<3> checker;
		checker.checkBasisFunctions(&basis_func, evaluation_points);

		// Derivatives
		std::vector<VectorDouble> derivatives;
		derivatives = basis_func.ComputeBasisFunctionDerivatives(onezerozero);
		TS_ASSERT_EQUALS(derivatives.size(), 10);
		TS_ASSERT_DELTA(derivatives[0](0),  1, 1e-12);
		TS_ASSERT_DELTA(derivatives[1](0),  3, 1e-12);
		TS_ASSERT_DELTA(derivatives[2](0),  0, 1e-12);
		TS_ASSERT_DELTA(derivatives[3](0),  0, 1e-12);
	}
	
	void TestComputeTransformedQuadraticBasisFunctionDerivatives1d( void )
	{
		std::vector<const Node<1>*> nodes;
		nodes.push_back(new Node<1>(0, false, 3.0));
		nodes.push_back(new Node<1>(1, false, 5.0));
		nodes.push_back(new Node<1>(2, false, 4.0));
		Element<1,1> element(nodes,2);
		QuadraticBasisFunction<1> basis_function;
		
		const MatrixDouble *inverseJacobian = element.GetInverseJacobian();
		Point<1> evaluation_point(0.2); 
		std::vector<VectorDouble> trans_deriv =
			basis_function.ComputeTransformedBasisFunctionDerivatives(evaluation_point,
																		*inverseJacobian);
																		
		TS_ASSERT_DELTA(trans_deriv[0](0),-1.1, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[1](0),-0.1, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[2](0),1.2, 1e-12);
	}
	
	void TestComputeTransformedQuadraticBasisFunction2d( void )		
	{		
		std::vector<const Node<2>*> nodes;
		nodes.push_back(new Node<2>(0, false, 4.0, 3.0));
		nodes.push_back(new Node<2>(1, false, 6.0, 4.0));
		nodes.push_back(new Node<2>(2, false, 3.0, 5.0));
		nodes.push_back(new Node<2>(3, false, 5.0, 3.5));
		nodes.push_back(new Node<2>(4, false, 3.5, 4.0));
		nodes.push_back(new Node<2>(5, false, 4.5, 4.5));
		Element<2,2> element(nodes,2);
		QuadraticBasisFunction<2> basis_function;
		
		const MatrixDouble *inverseJacobian = element.GetInverseJacobian();
		Point<2> evaluation_point(0.3, 0.6); 
		std::vector<VectorDouble> trans_deriv =
			basis_function.ComputeTransformedBasisFunctionDerivatives(evaluation_point,
																		*inverseJacobian);
		
		TS_ASSERT_DELTA(trans_deriv[0](0),0.12, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[0](1),0.36, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[1](0),0.08, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[1](1),0.04, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[2](0),-0.28, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[2](1),0.56, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[3](0),-0.08, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[3](1),-0.64, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[4](0),-0.56, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[4](1),-1.28, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[5](0),0.72, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[5](1),0.96, 1e-12);
	}
	
	void TestComputeTransformedQuadraticBasisFunction3d( void )		
	{		
		std::vector<const Node<3>*> nodes;
		nodes.push_back(new Node<3>(0, false, 4.0, 3.0, 0.0));
		nodes.push_back(new Node<3>(1, false, 6.0, 4.0, 1.0));
		nodes.push_back(new Node<3>(2, false, 3.0, 5.0, 2.0));
		nodes.push_back(new Node<3>(3, false, 5.0, 4.0, 3.0));
		nodes.push_back(new Node<3>(4, false, 5.0, 3.5, 0.5));
		nodes.push_back(new Node<3>(5, false, 3.5, 4.0, 1.0));
		nodes.push_back(new Node<3>(6, false, 4.5, 3.5, 1.5));
		nodes.push_back(new Node<3>(7, false, 4.5, 4.5, 1.5));
		nodes.push_back(new Node<3>(8, false, 5.5, 4.0, 2.0));
		nodes.push_back(new Node<3>(9, false, 4.0, 4.5, 2.5));
		Element<3,3> element(nodes,2);
		QuadraticBasisFunction<3> basis_function;
		
		const MatrixDouble *inverseJacobian = element.GetInverseJacobian();
		Point<3> evaluation_point(0.3, 0.1, 0.2); 
		std::vector<VectorDouble> trans_deriv =
			basis_function.ComputeTransformedBasisFunctionDerivatives(evaluation_point,
																		*inverseJacobian);
		
		TS_ASSERT_DELTA(trans_deriv[0](0),-0.12, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[0](1),-0.3, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[0](2),-0.06, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[1](0),0.08, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[1](1),0.1, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[1](2),-0.06, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[2](0),0.12, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[2](1),-0.3, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[2](2),0.06, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[3](0),0.0, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[3](1),0.1, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[3](2),-0.1, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[4](0),0.4, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[4](1),0.2, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[4](2),-0.6, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[5](0),-0.4, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[5](1),0.6, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[5](2),-0.2, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[6](0),-0.16, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[6](1),-1.2, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[6](2),0.72, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[7](0),-0.08, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[7](1),0.8, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[7](2),-0.24, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[8](0),0.32, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[8](1),-0.2, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[8](2),0.36, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[9](0),-0.16, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[9](1),0.2, 1e-12);
		TS_ASSERT_DELTA(trans_deriv[9](2),0.12, 1e-12);
	}
};

#endif //_TESTQUADRATICBASISFUNCTION_HPP_
