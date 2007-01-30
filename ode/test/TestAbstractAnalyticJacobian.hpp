#ifndef _TESTABSTRACTANALYTICJACOBIAN_HPP_
#define _TESTABSTRACTANALYTICJACOBIAN_HPP_

// TestAbstractAnalyticJacobian.hpp

#include <cmath>
#include <iostream>
#include <vector>
#include "Jacobian1.hpp"
#include "Jacobian2.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PetscException.hpp"

// Tolerance for tests
double tol=0.01;


class TestAbstractAnalyticJacobian : public CxxTest::TestSuite
{
public:

    void TestJacobianOne(void)
    {
        // pointer to TestOde1 class
        Jacobian1 ode_system;
        
        std::vector<double>  solution_guess(1);
        solution_guess[0] = 2.0;
        
        // Set up a Jacobian matrix for function to put values in
        double** jacobian;
        
        jacobian = new double*[1];
        jacobian[0] = new double[1];
        
        // This is the function we are testing...
        ode_system.BetterAnalyticJacobian(solution_guess, jacobian, 1.0, 0.01);
        
        TS_ASSERT_DELTA(jacobian[0][0], 0.96, tol);
             
        delete jacobian[0];
        delete jacobian; 
    }
    
    void TestJacobianTwo(void)
    {
        // pointer to TestOde1 class
        Jacobian2 ode_system;
        
        std::vector<double>  solution_guess(2);
        solution_guess[0] = 1.0;
        solution_guess[1] = 2.0;
        
        // Set up a Jacobian matrix for function to put values in
        double** jacobian;
        
        jacobian = new double*[2];
        jacobian[0] = new double[2];
        jacobian[1] = new double[2];
        
        // This is the function we are testing...
        ode_system.BetterAnalyticJacobian(solution_guess, jacobian, 1.0, 0.01);
        
        TS_ASSERT_DELTA(jacobian[0][0], 0.98, tol);
        TS_ASSERT_DELTA(jacobian[0][1], -0.04, tol);
        TS_ASSERT_DELTA(jacobian[1][0], -0.02, tol);
        TS_ASSERT_DELTA(jacobian[1][1], 0.92, tol);             
    }
        
};



#endif //_TESTABSTRACTANALYTICJACOBIAN_HPP_
