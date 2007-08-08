#ifndef _TESTLINEARBASISFUNCTION_HPP_
#define _TESTLINEARBASISFUNCTION_HPP_

#include <cxxtest/TestSuite.h>
#include <vector>

#include "LinearBasisFunction.cpp"
#include "GaussianQuadratureRule.hpp"
#include "BasisFunctionsCheckers.hpp"
#include "Element.hpp"

class TestLinearBasisFunction : public CxxTest::TestSuite
{
public:

    void TestLinearBasisFunction0d()
    {
        ChastePoint<0> zero;
        TS_ASSERT_DELTA(LinearBasisFunction<0>::ComputeBasisFunction(zero, 0), 1.0, 1e-12);
        
        c_vector<double, 1> basis_function_vector;
        basis_function_vector = LinearBasisFunction<0>::ComputeBasisFunctions(zero);
        TS_ASSERT_EQUALS(basis_function_vector.size(), 1u);
        TS_ASSERT_DELTA(basis_function_vector[0], 1.0, 1e-12);
        
        // check link with 0d quad rule works ok
        GaussianQuadratureRule<0>  quad_rule(1);
        const ChastePoint<0>& quad_point = quad_rule.rGetQuadPoint(0);
        
        c_vector<double, 1> basis_function_vector2;
        basis_function_vector2 = LinearBasisFunction<0>::ComputeBasisFunctions(quad_point);
        TS_ASSERT_EQUALS(basis_function_vector.size(), 1u);
        TS_ASSERT_DELTA(basis_function_vector[0], 1.0, 1e-12);
    }
    
    void TestLinearBasisFunction1d()
    {
        ChastePoint<1> zero(0);
        ChastePoint<1> one(1);
        
        std::vector<ChastePoint<1>*> evaluation_points;
        evaluation_points.push_back(&zero);
        evaluation_points.push_back(&one);
        
        BasisFunctionsCheckers<1> checker;
        checker.checkBasisFunctions(evaluation_points);
        
        // Derivatives
        c_matrix<double, 1, 2> derivatives;
        derivatives = LinearBasisFunction<1>::ComputeBasisFunctionDerivatives(one);
        TS_ASSERT_DELTA(derivatives(0,0), -1, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,1),  1, 1e-12);
    }
    
    void TestLinearBasisFunction2d()
    {
        ChastePoint<2> zero(0,0);
        ChastePoint<2> onezero(1,0);
        ChastePoint<2> zeroone(0,1);
        
        std::vector<ChastePoint<2>*> evaluation_points;
        evaluation_points.push_back(&zero);
        evaluation_points.push_back(&onezero);
        evaluation_points.push_back(&zeroone);
        
        BasisFunctionsCheckers<2> checker;
        checker.checkBasisFunctions(evaluation_points);
        
        // Derivatives
        c_matrix<double, 2, 3> derivatives;
        derivatives = LinearBasisFunction<2>::ComputeBasisFunctionDerivatives(onezero);
        TS_ASSERT_DELTA(derivatives(0,0), -1, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,1),  1, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,2),  0, 1e-12);
    }
    
    
    void TestLinearBasisFunction3d()
    {
        ChastePoint<3> zero(0,0,0);
        ChastePoint<3> zerozeroone(0,0,1);
        ChastePoint<3> zeroonezero(0,1,0);
        ChastePoint<3> onezerozero(1,0,0);
        
        std::vector<ChastePoint<3>*> evaluation_points;
        evaluation_points.push_back(&zero);
        evaluation_points.push_back(&onezerozero);
        evaluation_points.push_back(&zeroonezero);
        evaluation_points.push_back(&zerozeroone);
        
        BasisFunctionsCheckers<3> checker;
        checker.checkBasisFunctions(evaluation_points);
        
        // Derivatives
        c_matrix<double, 3, 4> derivatives;
        derivatives = LinearBasisFunction<3>::ComputeBasisFunctionDerivatives(onezerozero);
        TS_ASSERT_DELTA(derivatives(0,0), -1, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,1),  1, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,2),  0, 1e-12);
        TS_ASSERT_DELTA(derivatives(0,3),  0, 1e-12);
    }
    
    
    void TestComputeTransformedBasisFunctionDerivatives( void )
    {
        // 1D
        ChastePoint<1> one(1);
        
        c_matrix<double, 1, 1> inv_J;
        inv_J(0,0)=0.5;
        
        c_matrix<double, 1, 2> transDeriv =
            LinearBasisFunction<1>::ComputeTransformedBasisFunctionDerivatives(one, inv_J);
            
        TS_ASSERT_DELTA(transDeriv(0,0), -0.5, 1e-12);
        TS_ASSERT_DELTA(transDeriv(0,1),  0.5, 1e-12);
        
        // 2D
        ChastePoint<2> oneone(1,1);
        
        c_matrix<double, 2, 2> inv_J2 = 0.5 * identity_matrix<double>(2);
        
        c_matrix<double, 2, 3> trans_deriv =
            LinearBasisFunction<2>::ComputeTransformedBasisFunctionDerivatives(oneone, inv_J2);
            
        TS_ASSERT_DELTA(trans_deriv(0,0), -0.5, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(0,1),  0.5, 1e-12);
        TS_ASSERT_DELTA(trans_deriv(0,2),    0, 1e-12);
        
        //3D
        ChastePoint<3> oneoneone(1,1,1);
        
        c_matrix<double, 3, 3> inv_J3 = 0.5 * identity_matrix<double>(3);
        
        c_matrix<double, 3, 4> transDeriv3 =
            LinearBasisFunction<3>::ComputeTransformedBasisFunctionDerivatives(oneoneone,inv_J3);
            
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
        
        const c_matrix<double, 2, 2> *inverseJacobian = element.GetInverseJacobian();
        ChastePoint<2> evaluation_point(1,1);
        c_matrix<double, 2, 3> trans_deriv =
            LinearBasisFunction<2>::ComputeTransformedBasisFunctionDerivatives(evaluation_point,
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
