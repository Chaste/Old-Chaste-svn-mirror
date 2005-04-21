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
	
};

#endif //_TESTLINEARBASISFUNCTION_HPP_
