#ifndef _TESTLINEARBASISFUNCTION_HPP_
#define _TESTLINEARBASISFUNCTION_HPP_

#include <cxxtest/TestSuite.h>
#include "LinearBasisFunction.cpp"
#include "GaussianQuadratureRule.hpp"
#include "BasisFunctionsCheckers.hpp"
#include "Element.hpp"
#include <vector>

class TestLinearBasisFunction : public CxxTest::TestSuite
{
public:

    void TestLinearBasisFunction0d()
    {
        Point<0> zero;
        LinearBasisFunction<0> basis_function;
        TS_ASSERT_DELTA(basis_function.ComputeBasisFunction(zero, 0), 1.0, 1e-12);
        
        c_vector<double, 1> basis_function_vector;
        basis_function_vector = basis_function.ComputeBasisFunctions(zero);
        TS_ASSERT_EQUALS(basis_function_vector.size(), 1u);
        TS_ASSERT_DELTA(basis_function_vector[0], 1.0, 1e-12);
        
        // check link with 0d quad rule works ok
        GaussianQuadratureRule<0>  quad_rule(1);
        Point<0>   quad_point = quad_rule.GetQuadPoint(0);
        
        c_vector<double, 1> basis_function_vector2;
        basis_function_vector2 = basis_function.ComputeBasisFunctions(quad_point);
        TS_ASSERT_EQUALS(basis_function_vector.size(), 1u);
        TS_ASSERT_DELTA(basis_function_vector[0], 1.0, 1e-12);
    }
    
    void TestLinearBasisFunction1d()
    {
        Point<1> zero(0);
        Point<1> one(1);
        LinearBasisFunction<1> basis_func;
        
        std::vector<Point<1>*> evaluation_points;
        evaluation_points.push_back(&zero);
        evaluation_points.push_back(&one);
        
        BasisFunctionsCheckers<1> checker;
        checker.checkBasisFunctions(&basis_func, evaluation_points);
        
        // Derivatives
        c_matrix<double, 1, 2> derivatives;
        derivatives = basis_func.ComputeBasisFunctionDerivatives(one);
        TS_ASSERT_DELTA(derivatives(0,0), -1, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,1),  1, 1e-12);
    }
    
    void TestLinearBasisFunction2d()
    {
        Point<2> zero(0,0);
        Point<2> onezero(1,0);
        Point<2> zeroone(0,1);
        
        LinearBasisFunction<2> basis_func;
        
        std::vector<Point<2>*> evaluation_points;
        evaluation_points.push_back(&zero);
        evaluation_points.push_back(&onezero);
        evaluation_points.push_back(&zeroone);
        
        BasisFunctionsCheckers<2> checker;
        checker.checkBasisFunctions(&basis_func, evaluation_points);
        
        // Derivatives
        c_matrix<double, 2, 3> derivatives;
        derivatives = basis_func.ComputeBasisFunctionDerivatives(onezero);
        TS_ASSERT_DELTA(derivatives(0,0), -1, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,1),  1, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,2),  0, 1e-12);
    }
    
    
    void TestLinearBasisFunction3d()
    {
        Point<3> zero(0,0,0);
        Point<3> zerozeroone(0,0,1);
        Point<3> zeroonezero(0,1,0);
        Point<3> onezerozero(1,0,0);
        
        LinearBasisFunction<3> basis_func;
        
        std::vector<Point<3>*> evaluation_points;
        evaluation_points.push_back(&zero);
        evaluation_points.push_back(&onezerozero);
        evaluation_points.push_back(&zeroonezero);
        evaluation_points.push_back(&zerozeroone);
        
        BasisFunctionsCheckers<3> checker;
        checker.checkBasisFunctions(&basis_func, evaluation_points);
        
        // Derivatives
        c_matrix<double, 3, 4> derivatives;
        derivatives = basis_func.ComputeBasisFunctionDerivatives(onezerozero);
        TS_ASSERT_DELTA(derivatives(0,0), -1, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,1),  1, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,2),  0, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,3),  0, 1e-12);
    }
    
    
    void TestComputeTransformedBasisFunctionDerivatives( void )
    {
        // 1D
        LinearBasisFunction<1> basis_func;
        Point<1> one(1);
        
        c_matrix<double, 1, 1> inv_J;
        inv_J(0,0)=0.5;
        
        c_matrix<double, 1, 2> transDeriv =
            basis_func.ComputeTransformedBasisFunctionDerivatives(one, inv_J);
            
        TS_ASSERT_DELTA(transDeriv(0,0), -0.5, 1e-12);
        TS_ASSERT_DELTA(transDeriv(0,1),  0.5, 1e-12);
        
        // 2D
        LinearBasisFunction<2> basis_func2;
        Point<2> oneone(1,1);
        
        c_matrix<double, 2, 2> inv_J2 = 0.5 * identity_matrix<double>(2);
        
        c_matrix<double, 2, 3> trans_deriv =
            basis_func2.ComputeTransformedBasisFunctionDerivatives(oneone, inv_J2);
            
        TS_ASSERT_DELTA(trans_deriv(0,0), -0.5, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(0,1),  0.5, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(0,2),    0, 1e-12);
        
        
        //3D
        LinearBasisFunction<3> basis_func3;
        Point<3> oneoneone(1,1,1);
        
        c_matrix<double, 3, 3> inv_J3 = 0.5 * identity_matrix<double>(3);
        
        c_matrix<double, 3, 4> transDeriv3 =
            basis_func3.ComputeTransformedBasisFunctionDerivatives(oneoneone,inv_J3);
            
        TS_ASSERT_DELTA(transDeriv3(0,0), -0.5, 1e-12);
        TS_ASSERT_DELTA(transDeriv3(0,1),  0.5, 1e-12);
        TS_ASSERT_DELTA(transDeriv3(0,2),    0, 1e-12);
        TS_ASSERT_DELTA(transDeriv3(0,3),    0, 1e-12);
        //TS_TRACE("here lin basis\n");
    }
    
    void TestComputeTransformedBasisFunction2( void )
    {
        // 2D - with better test data
        
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 4.0, 3.0));
        nodes.push_back(new Node<2>(1, false, 6.0, 4.0));
        nodes.push_back(new Node<2>(2, false, 3.0, 5.0));
        Element<2,2> element(INDEX_IS_NOT_USED, nodes);
        LinearBasisFunction<2> basis_function;
        
        const c_matrix<double, 2, 2> *inverseJacobian = element.GetInverseJacobian();
        Point<2> evaluation_point(1,1);
        c_matrix<double, 2, 3> trans_deriv =
            basis_function.ComputeTransformedBasisFunctionDerivatives(evaluation_point,
                                                                      *inverseJacobian);
                                                                      
        TS_ASSERT_DELTA(trans_deriv(0,0),-0.2, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(1,0),-0.6, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(0,1),0.4, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(1,1),0.2, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(0,2),-0.2, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(1,2),0.4, 1e-12);
        
        delete nodes[0];
        delete nodes[1];
        delete nodes[2];        
    }
};

#endif //_TESTLINEARBASISFUNCTION_HPP_
