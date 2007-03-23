#ifndef _TESTABSTRACTODESYSTEM_HPP_
#define _TESTABSTRACTODESYSTEM_HPP_

// TestAbstractOdeSystem.hpp

#include <cmath>
#include <iostream>
#include <vector>
#include "AbstractOdeSystem.hpp"
#include "Ode1.hpp"
#include "Ode2.hpp"
#include "Ode3.hpp"
#include "TwoDimOdeSystem.hpp"

// Tolerance for tests
double tol=0.01;


class TestAbstractOdeSystem : public CxxTest::TestSuite
{
public:

    void TestOdeSystemOne(void)
    {
        // Test Ode1 class
        Ode1 ode1;
        // dy
        std::vector<double> dy(1);
        ode1.EvaluateYDerivatives(1.0, ode1.GetInitialConditions(),dy);
        TS_ASSERT_DELTA(dy[0],1.0,tol);
    }
    
    
    void TestOdeSystemTwo(void)
    {
        Ode2 ode2;
        std::vector<double> dy(1);
        ode2.EvaluateYDerivatives(2.0, ode2.GetInitialConditions(), dy);
        TS_ASSERT_DELTA(dy[0],8.0,tol);
    }
    
    
    void TestOdeSystemThree(void)
    {
        Ode3 ode3;
        std::vector<double> dy(2);
        ode3.EvaluateYDerivatives(2.0, ode3.GetInitialConditions(), dy);
        TS_ASSERT_DELTA(dy[0],8.0,tol);
        TS_ASSERT_DELTA(dy[1],16.0,tol);
    }


    void TestExceptions(void)
    {
        Ode1 ode;
        TS_ASSERT_EQUALS(ode.GetNumberOfStateVariables(), 1u);
        std::vector<double> v(2);
        v[0] = -1.0;
        v[1] = -2.0;
        TS_ASSERT_THROWS_ANYTHING(ode.SetInitialConditions(v));
        TS_ASSERT_THROWS_ANYTHING(ode.SetInitialConditionsComponent(2, -3.0));
        TS_ASSERT_THROWS_ANYTHING(ode.SetStateVariables(v));
    }
    
    
    void TestSetGetFunctionsInAbstractOdeSystem(void)
    {
        TwoDimOdeSystem ode;
        
        std::vector<double> initial_conditions = ode.GetInitialConditions();
        std::vector<double> state_variables = ode.rGetStateVariables();
        
        TS_ASSERT_DELTA( initial_conditions[0], 1.0, 1e-12 );
        TS_ASSERT_DELTA( initial_conditions[1], 2.0, 1e-12 );
        TS_ASSERT_DELTA( state_variables[0], 3.0, 1e-12 );
        TS_ASSERT_DELTA( state_variables[1], 4.0, 1e-12 );
        
        std::vector<double> new_initial_conditions;
        new_initial_conditions.push_back(5.0);
        new_initial_conditions.push_back(6.0);
        
        std::vector<double> new_state_variables;
        new_state_variables.push_back(7.0);
        new_state_variables.push_back(8.0);
        
        ode.SetInitialConditions(new_initial_conditions);
        ode.SetStateVariables(new_state_variables);
        
        initial_conditions = ode.GetInitialConditions();
        state_variables = ode.rGetStateVariables();
        
        TS_ASSERT_DELTA( initial_conditions[0], 5.0, 1e-12 );
        TS_ASSERT_DELTA( initial_conditions[1], 6.0, 1e-12 );
        TS_ASSERT_DELTA( state_variables[0], 7.0, 1e-12 );
        TS_ASSERT_DELTA( state_variables[1], 8.0, 1e-12 );
        
        ode.SetInitialConditionsComponent(1,9.0);
        initial_conditions = ode.GetInitialConditions();
        TS_ASSERT_DELTA( initial_conditions[0], 5.0, 1e-12 );
        TS_ASSERT_DELTA( initial_conditions[1], 9.0, 1e-12 );
    }
};



#endif //_TESTABSTRACTODESYSTEM_HPP_
