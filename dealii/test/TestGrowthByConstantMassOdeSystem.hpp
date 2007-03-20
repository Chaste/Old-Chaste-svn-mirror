#ifndef TESTGROWTHBYCONSTANTMASSODESYSTEM_HPP_
#define TESTGROWTHBYCONSTANTMASSODESYSTEM_HPP_

#include <cxxtest/TestSuite.h>

#include "GrowthByConstantMassOdeSystem.hpp"
#include "SimpleTumourSourceModel.hpp"
#include "EulerIvpOdeSolver.hpp"

class TestGrowthByConstantMassOdeSystem : public CxxTest::TestSuite
{
public:
    void TestWithSimpleTumourSourceModel2d() throw(Exception)
    {
        SimpleTumourSourceModel<2> source_model;
        
        Point<2> position;
        double rho = 1.0;
        
        source_model.AddEvaluationPoint(0,position,0);
        GrowthByConstantMassOdeSystem<2> ode_system_0(rho, 0, &source_model);
        
        source_model.AddEvaluationPoint(1,position,0);
        GrowthByConstantMassOdeSystem<2> ode_system_1(rho, 1, &source_model);
        
        double start_time = 0;
        double end_time = 1;
        double dt = 0.01;
        source_model.Run(start_time, end_time, NULL);
        
        EulerIvpOdeSolver solver;
        solver.SolveAndUpdateStateVariable(&ode_system_0, start_time, end_time, dt);
        
        //rho=1, s=0, so ode is dg/dt = 0, soln is g=const
        //initial condition g=1
        TS_ASSERT_DELTA(ode_system_0.rGetStateVariables()[0], 1.0, 1e-2);
        
        solver.SolveAndUpdateStateVariable(&ode_system_1, start_time, end_time, dt);
        
        //rho=1, s=1, so ode is dg/dt = (1/2)g, soln is g=Ae^{0.5t}
        //initial condition => A=1, so at t=1, g = e^{0.5} = 1.64872127
        double sqrt_e = 1.64872127;
        TS_ASSERT_DELTA(ode_system_1.rGetStateVariables()[0], sqrt_e, 1e-2);
    }
    
    
    
    void TestWithSimpleTumourSourceModel3d() throw(Exception)
    {
        SimpleTumourSourceModel<3> source_model;
        
        Point<3> position;
        double rho = 1.0;
        
        unsigned eval_point_index = 1;
        source_model.AddEvaluationPoint(eval_point_index,position,0);
        GrowthByConstantMassOdeSystem<3> ode_system(rho, eval_point_index, &source_model);
        
        double start_time = 0;
        double end_time = 1;
        double dt = 0.01;
        source_model.Run(start_time, end_time, NULL);
        
        EulerIvpOdeSolver solver;
        solver.SolveAndUpdateStateVariable(&ode_system, start_time, end_time, dt);
        
        //rho=1, s=1, so ode is dg/dt = (1/3)g, soln is g=Ae^{t/3}
        //initial condition => A=1, so at t=1, g = e^{1/3} = 1.39561243
        double e_to_the_third = 1.39561243;
        TS_ASSERT_DELTA(ode_system.rGetStateVariables()[0], e_to_the_third, 1e-2);
    }
};



#endif /*TESTGROWTHBYCONSTANTMASSODESYSTEM_HPP_*/
