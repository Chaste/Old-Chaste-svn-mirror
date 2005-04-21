#ifndef _TESTLINEARBASISFUNCTION_HPP_
#define _TESTLINEARBASISFUNCTION_HPP_

#include <cxxtest/TestSuite.h>
#include "LinearBasisFunction.cpp"
#include "VectorDouble.hpp"
#include <vector>

class TestLinearBasisFunction : public CxxTest::TestSuite 
{
	public:
	
	void testLinearBasisFunction1d()
	{
		Point<1> zero(0);
		Point<1> one(1);
		LinearBasisFunction<1> basis_func;
		
		// Single compute
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 0), 1.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(zero, 1), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(one, 0), 0.0, 1e-12);
		TS_ASSERT_DELTA(basis_func.ComputeBasisFunction(one, 1), 1.0, 1e-12);
		
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
		TS_ASSERT_DELTA(derivatives[1](0), 1, 1e-12);
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
		TS_ASSERT_DELTA(derivatives[1](0), 1, 1e-12);
		TS_ASSERT_DELTA(derivatives[2](0), 0, 1e-12);
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
		
		std::cout << " Nick Trefethen is Distressed !!  \n";
		// Derivatives
		std::vector<VectorDouble> derivatives;
		derivatives = basis_func.ComputeBasisFunctionDerivatives(onezerozero);
		TS_ASSERT_EQUALS(derivatives.size(), 4);
		TS_ASSERT_DELTA(derivatives[0](0), -1, 1e-12);
		TS_ASSERT_DELTA(derivatives[1](0), 1, 1e-12);
		TS_ASSERT_DELTA(derivatives[2](0), 0, 1e-12);
		TS_ASSERT_DELTA(derivatives[3](0), 0, 1e-12);
	}
	
};

#endif //_TESTLINEARBASISFUNCTION_HPP_
