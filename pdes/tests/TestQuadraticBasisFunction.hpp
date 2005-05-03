#ifndef _TESTQUADRATICBASISFUNCTION_HPP_
#define _TESTQUADRATICBASISFUNCTION_HPP_

#include <cxxtest/TestSuite.h>
#include "QuadraticBasisFunction.cpp"
#include "GaussianQuadratureRule.hpp"
#include "VectorDouble.hpp"
#include <vector>

class TestQuadraticBasisFunction : public CxxTest::TestSuite 
{
	public:
	
	void testQuadraticBasisFunction0d() 
	{
		Point<0> zero;
		QuadraticBasisFunction<0> basis_func;
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 0), 1.0, 1e-12);
	}
	
	void checkBasisFunctions1d(AbstractBasisFunction<1>* pBasisFunc,
							   std::vector<Point<1>*> evaluationPoints)
	{
		int size = evaluationPoints.size();		// number of evalutation points and basis functions too
		std::vector<double> basis_function_vector; // store results of evalutation
		double expected_evaluation;
		
		for (int point_index=0; point_index<size; point_index++)
		{
			std::vector<double> basis_function_vector = pBasisFunc->ComputeBasisFunctions(*(evaluationPoints[point_index]));
			for (int func_index=0; func_index<size; func_index ++)
			{
				if (func_index==point_index)
				{
					expected_evaluation=1.0;
				}
				else
				{
					expected_evaluation = 0.0;
				}
				TS_ASSERT_DELTA(basis_function_vector[func_index],
								expected_evaluation,
								1e-12);
				TS_ASSERT_DELTA(pBasisFunc->ComputeBasisFunction(*(evaluationPoints[point_index]),func_index),
								expected_evaluation,
								1e12);
			}
		}
	}
		
	
	void testQuadraticBasisFunction1d()
	{
		std::vector<Point<1>*> evaluation_points;
		Point<1> zero(0);
		Point<1> one(1);
		Point<1> half(0.5);
		evaluation_points.push_back(new Point<1>(0));
		evaluation_points.push_back(new Point<1>(0.5));
		evaluation_points.push_back(new Point<1>(1));
		
		QuadraticBasisFunction<1> basis_func;
		
		// checkBasisFunctions1d(&basis_func, evaluation_points);
		
		// Single compute
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 0), 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 1), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 2), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(half,  0), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(half,  1), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(half,  2), 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(one,  0), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(one,  1), 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(one,  2), 0.0, 1e-12);
		
		// Mass compute
		std::vector<double> basis_function_vector;
		basis_function_vector = basis_func.ComputeBasisFunctions(zero);
		TS_ASSERT_EQUALS(basis_function_vector.size(), 3);
		
		TS_ASSERT_DELTA(basis_function_vector[0], 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_function_vector[1], 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_function_vector[2], 0.0, 1e-12);
		
		basis_function_vector = basis_func.ComputeBasisFunctions(one);
		TS_ASSERT_DELTA(basis_function_vector[0], 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_function_vector[1], 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_function_vector[2], 0.0, 1e-12);
		
		basis_function_vector = basis_func.ComputeBasisFunctions(half);
		TS_ASSERT_DELTA(basis_function_vector[0], 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_function_vector[1], 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_function_vector[2], 1.0, 1e-12);
		
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
	
	void testQuadraticBasisFunction2d()
	{
		Point<2> zero(0,0);
		Point<2> onezero(1,0);
		Point<2> zeroone(0,1);
		Point<2> halfzero(0.5,0);
		Point<2> zerohalf(0,0.5);
		Point<2> halfhalf(0.5,0.5);
		
		QuadraticBasisFunction<2> basis_func;
		
		// Single compute
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 0), 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 1), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 2), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 3), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 4), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 5), 0.0, 1e-12);
		
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(onezero, 0), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(onezero, 1), 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(onezero, 2), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(onezero, 3), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(onezero, 4), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(onezero, 5), 0.0, 1e-12);
		
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zeroone, 0), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zeroone, 1), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zeroone, 2), 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zeroone, 3), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zeroone, 4), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zeroone, 5), 0.0, 1e-12);
		
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfzero, 0), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfzero, 1), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfzero, 2), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfzero, 3), 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfzero, 4), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfzero, 5), 0.0, 1e-12);
		
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerohalf, 0), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerohalf, 1), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerohalf, 2), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerohalf, 3), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerohalf, 4), 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerohalf, 5), 0.0, 1e-12);
		
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfhalf, 0), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfhalf, 1), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfhalf, 2), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfhalf, 3), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfhalf, 4), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfhalf, 5), 1.0, 1e-12);
		
		// Mass compute
		std::vector<double> basis_function_vector;
		basis_function_vector = basis_func.ComputeBasisFunctions(zero);
		TS_ASSERT_EQUALS(basis_function_vector.size(), 6);
		TS_ASSERT_DELTA(basis_function_vector[0], 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_function_vector[1], 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_function_vector[2], 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_function_vector[3], 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_function_vector[4], 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_function_vector[5], 0.0, 1e-12);
		
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


	void testQuadraticBasisFunction3d()
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
		
		QuadraticBasisFunction<3> basis_func;
		
		// Single compute
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 0), 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 1), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 2), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 3), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 4), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 5), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 6), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 7), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 8), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 9), 0.0, 1e-12);
		
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(onezerozero, 0), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(onezerozero, 1), 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(onezerozero, 2), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(onezerozero, 3), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(onezerozero, 4), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(onezerozero, 5), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(onezerozero, 6), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(onezerozero, 7), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(onezerozero, 8), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(onezerozero, 9), 0.0, 1e-12);
		
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zeroonezero, 0), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zeroonezero, 1), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zeroonezero, 2), 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zeroonezero, 3), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zeroonezero, 4), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zeroonezero, 5), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zeroonezero, 6), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zeroonezero, 7), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zeroonezero, 8), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zeroonezero, 9), 0.0, 1e-12);

		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerozeroone, 0), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerozeroone, 1), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerozeroone, 2), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerozeroone, 3), 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerozeroone, 4), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerozeroone, 5), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerozeroone, 6), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerozeroone, 7), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerozeroone, 8), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerozeroone, 9), 0.0, 1e-12);
		
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfhalfzero, 0), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfhalfzero, 1), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfhalfzero, 2), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfhalfzero, 3), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfhalfzero, 4), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfhalfzero, 5), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfhalfzero, 6), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfhalfzero, 7), 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfhalfzero, 8), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfhalfzero, 9), 0.0, 1e-12);
		
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfzerozero, 0), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfzerozero, 1), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfzerozero, 2), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfzerozero, 3), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfzerozero, 4), 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfzerozero, 5), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfzerozero, 6), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfzerozero, 7), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfzerozero, 8), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfzerozero, 9), 0.0, 1e-12);
		
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfzerohalf, 0), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfzerohalf, 1), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfzerohalf, 2), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfzerohalf, 3), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfzerohalf, 4), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfzerohalf, 5), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfzerohalf, 6), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfzerohalf, 7), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfzerohalf, 8), 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(halfzerohalf, 9), 0.0, 1e-12);
		
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerohalfhalf, 0), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerohalfhalf, 1), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerohalfhalf, 2), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerohalfhalf, 3), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerohalfhalf, 4), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerohalfhalf, 5), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerohalfhalf, 6), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerohalfhalf, 7), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerohalfhalf, 8), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerohalfhalf, 9), 1.0, 1e-12);
		
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerohalfzero, 0), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerohalfzero, 1), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerohalfzero, 2), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerohalfzero, 3), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerohalfzero, 4), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerohalfzero, 5), 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerohalfzero, 6), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerohalfzero, 7), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerohalfzero, 8), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerohalfzero, 9), 0.0, 1e-12);
		
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerozerohalf, 0), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerozerohalf, 1), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerozerohalf, 2), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerozerohalf, 3), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerozerohalf, 4), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerozerohalf, 5), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerozerohalf, 6), 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerozerohalf, 7), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerozerohalf, 8), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zerozerohalf, 9), 0.0, 1e-12);		
		
		// Mass compute
		std::vector<double> basis_function_vector;
		basis_function_vector = basis_func.ComputeBasisFunctions(zero);
		TS_ASSERT_EQUALS(basis_function_vector.size(), 10);
		TS_ASSERT_DELTA(basis_function_vector[0], 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_function_vector[1], 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_function_vector[2], 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_function_vector[3], 0.0, 1e-12);
		
		// Derivatives
		std::vector<VectorDouble> derivatives;
		derivatives = basis_func.ComputeBasisFunctionDerivatives(onezerozero);
		TS_ASSERT_EQUALS(derivatives.size(), 10);
		TS_ASSERT_DELTA(derivatives[0](0),  1, 1e-12);
		TS_ASSERT_DELTA(derivatives[1](0),  3, 1e-12);
		TS_ASSERT_DELTA(derivatives[2](0),  0, 1e-12);
		TS_ASSERT_DELTA(derivatives[3](0),  0, 1e-12);
	}
	
	void testComputeTransformedQuadraticBasisFunctionDerivatives1d( void )
	{
		std::vector<const Node<1>*> nodes;
		nodes.push_back(new Node<1>(0, false, 3.0));
		nodes.push_back(new Node<1>(1, false, 5.0));
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
	
	void testComputeTransformedQuadraticBasisFunction2d( void )		
	{		
		std::vector<const Node<2>*> nodes;
		nodes.push_back(new Node<2>(0, false, 4.0, 3.0));
		nodes.push_back(new Node<2>(1, false, 6.0, 4.0));
		nodes.push_back(new Node<2>(2, false, 3.0, 5.0));
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
	
	void testComputeTransformedQuadraticBasisFunction3d( void )		
	{		
		std::vector<const Node<3>*> nodes;
		nodes.push_back(new Node<3>(0, false, 4.0, 3.0, 0.0));
		nodes.push_back(new Node<3>(1, false, 6.0, 4.0, 1.0));
		nodes.push_back(new Node<3>(2, false, 3.0, 5.0, 2.0));
		nodes.push_back(new Node<3>(3, false, 5.0, 4.0, 3.0));
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
