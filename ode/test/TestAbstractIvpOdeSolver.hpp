#ifndef _TESTABSTRACTIVPODESOLVER_HPP_
#define _TESTABSTRACTIVPODESOLVER_HPP_
#include <cxxtest/TestSuite.h>

#include <vector>
#include <iostream>

#include "AbstractIvpOdeSolver.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "AdamsBashforthIvpOdeSolver.hpp"
#include "RungeKutta2IvpOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "Ode1.hpp"
#include "Ode4.hpp"
#include "OdeFirstOrder.hpp"
#include "OdeSecondOrder.hpp"
#include "OdeSecondOrderWithEvents.hpp"
#include "OdeThirdOrder.hpp"

#include "PetscSetupAndFinalize.hpp"

#define PI 3.14159265


class TestAbstractIvpOdeSolver: public CxxTest::TestSuite
{
private :
    void testGenericSolver(AbstractIvpOdeSolver& rSolver, double startTime, double endTime, double dt, double samplingTime)
    {
        // Initialise the instances of our ode system and solution classes
        Ode1 ode_system;
        OdeSolution solutions;
        
        // Solving the ode problem. Note that dt and the sampling time
        // are different
        std::vector<double> state_variables = ode_system.GetInitialConditions();
        solutions = rSolver.Solve(&ode_system, state_variables, startTime, endTime, dt, samplingTime);
        
        int num_timesteps = solutions.GetNumberOfTimeSteps();
        
        // the number of timesteps should be (just about) equal to
        // end_time/sampling_time = 2/0.01 = 200
        TS_ASSERT_DELTA(num_timesteps, (endTime-startTime)/samplingTime, 1);
        // also check the size of the data is correct
        TS_ASSERT_EQUALS(solutions.rGetSolutions().size(), (unsigned) (num_timesteps+1));
        
        int last = num_timesteps;
        
        // Test to solution is correct
        double testvalue = solutions.rGetSolutions()[last][0];
        
        // exact solution of Ode1 is y=t-t0
        TS_ASSERT_DELTA(testvalue, endTime-startTime, 0.01);
        
        // Test second version of Solve
        ode_system.SetStateVariables(ode_system.GetInitialConditions());
        state_variables = ode_system.rGetStateVariables();
        rSolver.Solve(&ode_system, state_variables, startTime, endTime, dt);
        TS_ASSERT_DELTA(state_variables[0], endTime-startTime, 0.01);
        
        // no stopping event was specified in the ODE, so check the
        // solver correctly states it didn't stop due to a
        // stopping event.
        TS_ASSERT_EQUALS(rSolver.StoppingEventOccured(), false);
    }
    
    
    // test a given solver on an ode which has a stopping event defined
    void testSolverOnOdesWithEvents(AbstractIvpOdeSolver& rSolver)
    {
        // ode which has solution y0 = cos(t), and stopping event y0<0,
        // ie should stop when t = pi/2;
        OdeSecondOrderWithEvents ode_with_events;
        
        OdeSolution solutions;
        std::vector<double> state_variables = ode_with_events.GetInitialConditions();
        solutions = rSolver.Solve(&ode_with_events, state_variables, 0.0, 2.0, 0.001, 0.001);
        
        int num_timesteps = solutions.GetNumberOfTimeSteps();
        
        // final time should be around pi/2
        TS_ASSERT_DELTA( solutions.rGetTimes()[num_timesteps], PI/2, 0.01);
        
        // penultimate y0 should be greater than zero
        TS_ASSERT_LESS_THAN( 0, solutions.rGetSolutions()[num_timesteps-1][0]);
        
        // final y0 should be less than zero
        TS_ASSERT_LESS_THAN( solutions.rGetSolutions()[num_timesteps][0], 0);
        
        // solver should correctly state the stopping event occured
        TS_ASSERT_EQUALS(rSolver.StoppingEventOccured(), true);
        
        
        ///////////////////////////////////////////////
        // repeat with sampling time larger than dt
        ///////////////////////////////////////////////
        
        state_variables = ode_with_events.GetInitialConditions();
        solutions = rSolver.Solve(&ode_with_events, state_variables, 0.0, 2.0, 0.001, 0.01);
        
        num_timesteps = solutions.GetNumberOfTimeSteps();
        
        // final time should be around pi/2
        TS_ASSERT_DELTA( solutions.rGetTimes()[num_timesteps], PI/2, 0.01);
        
        // penultimate y0 should be greater than zero
        TS_ASSERT_LESS_THAN( 0, solutions.rGetSolutions()[num_timesteps-1][0]);
        
        // final y0 should be less than zero
        TS_ASSERT_LESS_THAN( solutions.rGetSolutions()[num_timesteps][0], 0);
        
        // solver should correctly state the stopping event occured
        TS_ASSERT_EQUALS(rSolver.StoppingEventOccured(), true);
        
        // cover the check event isn't initially true exception
        std::vector<double> bad_init_cond;
        bad_init_cond.push_back(-1);  //y0 < 0 so stopping event true
        bad_init_cond.push_back(0.0);
        TS_ASSERT_THROWS_ANYTHING(rSolver.Solve(&ode_with_events, bad_init_cond, 0.0, 2.0, 0.001, 0.01));
    }
    
    
    
    
    
public:

    void testEulerSolver()
    {
        EulerIvpOdeSolver euler_solver;
        
        testGenericSolver(euler_solver,  0.0, 2.0, 0.001, 0.001);
        testGenericSolver(euler_solver,  1.0, 2.0, 0.001, 0.01);
        testGenericSolver(euler_solver, -1.0, 2.0, 0.001, 2);
        testGenericSolver(euler_solver,  0.0, 0.4, 0.01,  0.34);
        
        testSolverOnOdesWithEvents(euler_solver);
    }
    
    void testAdamsBashforthSolver()
    {
        AdamsBashforthIvpOdeSolver adams_bashforth_solver;
        
        testGenericSolver(adams_bashforth_solver,  0.0, 2.0, 0.001, 0.001);
        testGenericSolver(adams_bashforth_solver,  1.0, 2.0, 0.001, 0.01);
        testGenericSolver(adams_bashforth_solver, -1.0, 2.0, 0.001, 2);
        testGenericSolver(adams_bashforth_solver,  0.0, 0.4, 0.01,  0.34);
        
        testSolverOnOdesWithEvents(adams_bashforth_solver);
        
        // check exception thrown if number of timesteps <= 4
        TS_ASSERT_THROWS_ANYTHING( testGenericSolver(adams_bashforth_solver,0.0,0.04,0.01,0.01) );
    }
    
    void testRungeKutta2Solver()
    {
        RungeKutta2IvpOdeSolver rk2_solver;
        
        testGenericSolver(rk2_solver,  0.0, 2.0, 0.001, 0.001);
        testGenericSolver(rk2_solver,  1.0, 2.0, 0.001, 0.01);
        testGenericSolver(rk2_solver, -1.0, 2.0, 0.001, 2);
        testGenericSolver(rk2_solver,  0.0, 0.4, 0.01,  0.34);
        
        testSolverOnOdesWithEvents(rk2_solver);
    }
    
    void testRungeKutta4Solver()
    {
        RungeKutta4IvpOdeSolver rk4_solver;
        
        testGenericSolver(rk4_solver,  0.0, 2.0, 0.001, 0.001);
        testGenericSolver(rk4_solver,  1.0, 2.0, 0.001, 0.01);
        testGenericSolver(rk4_solver, -1.0, 2.0, 0.001, 2);
        testGenericSolver(rk4_solver,  0.0, 0.4, 0.01,  0.34);
        
        testSolverOnOdesWithEvents(rk4_solver);
    }
    
    
    void testLastTimeStep()
    {
        Ode1 ode_system;
        
        // Initialise the instance of our solver class
        EulerIvpOdeSolver euler_solver;
        // Initialise the instance of our solution class
        OdeSolution solutions;
        
        // Solving the ode problem. Note that dt and the sampling time
        // are different
        std::vector<double> state_variables = ode_system.GetInitialConditions();
        solutions = euler_solver.Solve(&ode_system, state_variables, 0.0, 2.0, 0.000037, 0.000037);
        
        int last = solutions.GetNumberOfTimeSteps();
        // Test to see if this worked
        double testvalue = solutions.rGetSolutions()[last-1][0]	;
        
        TS_ASSERT_DELTA(testvalue,2.0,0.001);
    }
    
    
    
    void testGlobalError()
    {
        OdeFirstOrder ode_system;
        
        double h_value=0.01;
        
        //Euler solver solution worked out
        EulerIvpOdeSolver euler_solver;
        OdeSolution solutions_euler;
        
        std::vector<double> state_variables = ode_system.GetInitialConditions();
        solutions_euler = euler_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        int last = solutions_euler.GetNumberOfTimeSteps();
        double testvalue_euler = solutions_euler.rGetSolutions()[last][0];
        
        //Runge Kutta 2 solver solution worked out
        RungeKutta2IvpOdeSolver rk2_solver;
        OdeSolution solutions_rk2;
        
        state_variables = ode_system.GetInitialConditions();
        solutions_rk2 = rk2_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        int last2 = solutions_rk2.GetNumberOfTimeSteps();
        double testvalue_rk2 = solutions_rk2.rGetSolutions()[last2][0];
        
        //Runge Kutta 4 solver solution worked out
        RungeKutta4IvpOdeSolver rk4_solver;
        OdeSolution solutions_rk4;
        
        state_variables = ode_system.GetInitialConditions();
        solutions_rk4 = rk4_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        int last3 = solutions_rk4.GetNumberOfTimeSteps();
        double testvalue_rk4 = solutions_rk4.rGetSolutions()[last3][0];
        
        //Adams-Bashforth solver solution worked out
        AdamsBashforthIvpOdeSolver adams_bashforth_solver;
        OdeSolution solutions_adams_bashforth;
        
        state_variables = ode_system.GetInitialConditions();
        solutions_adams_bashforth = adams_bashforth_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        int last4 = solutions_adams_bashforth.GetNumberOfTimeSteps();
        double testvalue_adams_bashforth = solutions_adams_bashforth.rGetSolutions()[last4][0];
        
        // The tests
        double exact_solution=exp(2);
        
        double global_error_euler;
        global_error_euler = 0.5*exp(2)*(exp(2)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_euler,exact_solution,global_error_euler);
        
        double global_error_rk2;
        global_error_rk2 = (1.0/6.0)*h_value*exp(2)*(exp(2)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_rk2,exact_solution,global_error_rk2);
        
        double global_error_rk4;
        global_error_rk4 = (1.0/24.0)*pow(h_value,3)*exp(2)*(exp(2)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_rk4,exact_solution,global_error_rk4);
        
        double global_error_adams_bashforth;
        global_error_adams_bashforth = (1.0/6.0)*pow(h_value,3)*exp(2)*(exp(2)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_adams_bashforth,exact_solution,global_error_adams_bashforth);
    }
    
    
    void testGlobalErrorSystemOf2Equations()
    {
        OdeSecondOrder ode_system;
        
        double h_value=0.01;
        
        //Euler solver solution worked out
        EulerIvpOdeSolver euler_solver;
        OdeSolution solutions_euler;
        
        std::vector<double> state_variables = ode_system.GetInitialConditions();
        solutions_euler = euler_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        int last = solutions_euler.GetNumberOfTimeSteps();
        
        double testvalue_euler[2];
        testvalue_euler[0] = solutions_euler.rGetSolutions()[last][0];
        testvalue_euler[1] = solutions_euler.rGetSolutions()[last][1];
        
        //Runge Kutta 2 solver solution worked out
        RungeKutta2IvpOdeSolver rk2_solver;
        OdeSolution solutions_rk2;
        
        state_variables = ode_system.GetInitialConditions();
        solutions_rk2 = rk2_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        int last2 = solutions_rk2.GetNumberOfTimeSteps();
        
        double testvalue_rk2[2];
        testvalue_rk2[0] = solutions_rk2.rGetSolutions()[last2][0];
        testvalue_rk2[1] = solutions_rk2.rGetSolutions()[last2][1];
        
        //Runge Kutta 4 solver solution worked out
        RungeKutta4IvpOdeSolver rk4_solver;
        OdeSolution solutions_rk4;
        
        state_variables = ode_system.GetInitialConditions();
        solutions_rk4 = rk4_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        int last3 = solutions_rk4.GetNumberOfTimeSteps();
        double testvalue_rk4[2];
        testvalue_rk4[0] = solutions_rk4.rGetSolutions()[last3][0];
        testvalue_rk4[1] = solutions_rk4.rGetSolutions()[last3][1];
        
        //solutions_rk4.SaveToFile("result.dat");
        
        //Adams-Bashforth solver solution worked out
        AdamsBashforthIvpOdeSolver adams_bashforth_solver;
        OdeSolution solutions_adams_bashforth;
        
        state_variables = ode_system.GetInitialConditions();
        solutions_adams_bashforth = adams_bashforth_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        int last4 = solutions_adams_bashforth.GetNumberOfTimeSteps();
        double testvalue_adams_bashforth[2];
        testvalue_adams_bashforth[0] = solutions_adams_bashforth.rGetSolutions()[last4][0];
        testvalue_adams_bashforth[1] = solutions_adams_bashforth.rGetSolutions()[last4][1];
        
        // The tests
        double exact_solution[2];
        
        exact_solution[0] = sin(2);
        exact_solution[1] = cos(2);
        
        double global_error_euler;
        global_error_euler = 0.5*1*(exp(2)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_euler[0],exact_solution[0],global_error_euler);
        TS_ASSERT_DELTA(testvalue_euler[1],exact_solution[1],global_error_euler);
        
        double global_error_rk2;
        global_error_rk2 = (1.0/6.0)*h_value*1*(exp(2)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_rk2[0],exact_solution[0],global_error_rk2);
        TS_ASSERT_DELTA(testvalue_rk2[1],exact_solution[1],global_error_rk2);
        
        double global_error_rk4;
        global_error_rk4 = (1.0/24.0)*pow(h_value,3)*1*(exp(2)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_rk4[0],exact_solution[0],global_error_rk4);
        TS_ASSERT_DELTA(testvalue_rk4[1],exact_solution[1],global_error_rk4);
        
        double global_error_adams_bashforth;
        global_error_adams_bashforth = (1.0/6.0)*pow(h_value,3)*1*(exp(2)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_adams_bashforth[0],exact_solution[0],global_error_adams_bashforth);
        TS_ASSERT_DELTA(testvalue_adams_bashforth[1],exact_solution[1],global_error_adams_bashforth);
    }
    
    
    void testGlobalErrorSystemOf3Equations()
    {
        OdeThirdOrder ode_system;
        
        double h_value=0.01;
        
        //Euler solver solution worked out
        EulerIvpOdeSolver euler_solver;
        OdeSolution solutions_euler;
        
        std::vector<double> state_variables = ode_system.GetInitialConditions();
        solutions_euler = euler_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        int last = solutions_euler.GetNumberOfTimeSteps();
        
        double testvalue_euler[3];
        testvalue_euler[0] = solutions_euler.rGetSolutions()[last][0];
        testvalue_euler[1] = solutions_euler.rGetSolutions()[last][1];
        testvalue_euler[2] = solutions_euler.rGetSolutions()[last][2];
        
        //Runge Kutta 2 solver solution worked out
        RungeKutta2IvpOdeSolver rk2_solver;
        OdeSolution solutions_rk2;
        
        state_variables = ode_system.GetInitialConditions();
        solutions_rk2 = rk2_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        int last2 = solutions_rk2.GetNumberOfTimeSteps();
        
        double testvalue_rk2[3];
        testvalue_rk2[0] = solutions_rk2.rGetSolutions()[last2][0];
        testvalue_rk2[1] = solutions_rk2.rGetSolutions()[last2][1];
        testvalue_rk2[2] = solutions_rk2.rGetSolutions()[last2][2];
        
        //Runge Kutta 4 solver solution worked out
        RungeKutta4IvpOdeSolver rk4_solver;
        OdeSolution solutions_rk4;
        
        state_variables = ode_system.GetInitialConditions();
        solutions_rk4 = rk4_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        int last3 = solutions_rk4.GetNumberOfTimeSteps();
        
        double testvalue_rk4[3];
        testvalue_rk4[0] = solutions_rk4.rGetSolutions()[last3][0];
        testvalue_rk4[1] = solutions_rk4.rGetSolutions()[last3][1];
        testvalue_rk4[2] = solutions_rk4.rGetSolutions()[last3][2];
        
        //solutions_rk4.SaveToFile("result.dat");
        
        //Adams-Bashforth solver solution worked out
        AdamsBashforthIvpOdeSolver adams_bashforth_solver;
        OdeSolution solutions_adams_bashforth;
        
        state_variables = ode_system.GetInitialConditions();
        solutions_adams_bashforth = adams_bashforth_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        int last4 = solutions_adams_bashforth.GetNumberOfTimeSteps();
        
        double testvalue_adams_bashforth[3];
        testvalue_adams_bashforth[0] = solutions_adams_bashforth.rGetSolutions()[last4][0];
        testvalue_adams_bashforth[1] = solutions_adams_bashforth.rGetSolutions()[last4][1];
        testvalue_adams_bashforth[2] = solutions_adams_bashforth.rGetSolutions()[last4][2];
        
        // The tests
        double exact_solution[3];
        
        exact_solution[0] = -sin(2);
        exact_solution[1] = sin(2)+cos(2);
        exact_solution[2] = 2*sin(2);
        
        double global_error_euler;
        global_error_euler = 0.5*2*(exp(2)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_euler[0],exact_solution[0],global_error_euler);
        TS_ASSERT_DELTA(testvalue_euler[1],exact_solution[1],global_error_euler);
        TS_ASSERT_DELTA(testvalue_euler[2],exact_solution[2],global_error_euler);
        
        double global_error_rk2;
        global_error_rk2 = (1.0/6.0)*h_value*2*(exp(2)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_rk2[0],exact_solution[0],global_error_rk2);
        TS_ASSERT_DELTA(testvalue_rk2[1],exact_solution[1],global_error_rk2);
        TS_ASSERT_DELTA(testvalue_rk2[2],exact_solution[2],global_error_rk2);
        
        double global_error_rk4;
        global_error_rk4 = (1.0/24.0)*pow(h_value,3)*2*(exp(2)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_rk4[0],exact_solution[0],global_error_rk4);
        TS_ASSERT_DELTA(testvalue_rk4[1],exact_solution[1],global_error_rk4);
        TS_ASSERT_DELTA(testvalue_rk4[2],exact_solution[2],global_error_rk4);
        
        double global_error_adams_bashforth;
        global_error_adams_bashforth = (1.0/6.0)*pow(h_value,3)*2*(exp(2)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_adams_bashforth[0],exact_solution[0],global_error_adams_bashforth);
        TS_ASSERT_DELTA(testvalue_adams_bashforth[1],exact_solution[1],global_error_adams_bashforth);
        TS_ASSERT_DELTA(testvalue_adams_bashforth[2],exact_solution[2],global_error_adams_bashforth);
        
    }
    
    void testGlobalError2()
    {
        Ode4 ode_system;
        
        double h_value=0.001;
        
        //Euler solver solution worked out
        EulerIvpOdeSolver euler_solver;
        OdeSolution solutions_euler;
        
        std::vector<double> state_variables = ode_system.GetInitialConditions();
        solutions_euler = euler_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        int last = solutions_euler.GetNumberOfTimeSteps();
        double testvalue_euler = solutions_euler.rGetSolutions()[last][0];
        
        //Runge Kutta 2 solver solution worked out
        RungeKutta2IvpOdeSolver rk2_solver;
        OdeSolution solutions_rk2;
        
        state_variables = ode_system.GetInitialConditions();
        solutions_rk2 = rk2_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        int last2 = solutions_rk2.GetNumberOfTimeSteps();
        double testvalue_rk2 = solutions_rk2.rGetSolutions()[last2][0];
        
        //Runge Kutta 4 solver solution worked out
        RungeKutta4IvpOdeSolver rk4_solver;
        OdeSolution solutions_rk4;
        
        state_variables = ode_system.GetInitialConditions();
        solutions_rk4 = rk4_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        int last3 = solutions_rk4.GetNumberOfTimeSteps();
        double testvalue_rk4 = solutions_rk4.rGetSolutions()[last3][0];
        
        //Adams-Bashforth solver solution worked out
        AdamsBashforthIvpOdeSolver adams_bashforth_solver;
        OdeSolution solutions_adams_bashforth;
        
        state_variables = ode_system.GetInitialConditions();
        solutions_adams_bashforth = adams_bashforth_solver.Solve(&ode_system, state_variables, 0.0, 2.0, h_value, h_value);
        int last4 = solutions_adams_bashforth.GetNumberOfTimeSteps();
        double testvalue_adams_bashforth = solutions_adams_bashforth.rGetSolutions()[last4][0];
        
        // The tests
        double alpha = 100;
        double exact_solution=1/(1+exp(-alpha*2));
        
        double global_error_euler;
        global_error_euler = 0.5*1/(1+exp(-alpha*2))*(exp(2)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_euler,exact_solution,global_error_euler);
        
        double global_error_rk2;
        global_error_rk2 = (1.0/6.0)*h_value*1/(1+exp(-alpha*2))*(exp(2)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_rk2,exact_solution,global_error_rk2);
        
        double global_error_rk4;
        global_error_rk4 = (1.0/24.0)*pow(h_value,3)*1/(1+exp(-alpha*2))*(exp(2)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_rk4,exact_solution,global_error_rk4);
        
        double global_error_adams_bashforth;
        global_error_adams_bashforth = (1.0/6.0)*pow(h_value,3)*1/(1+exp(-alpha*2))*(exp(2)-1)*h_value;
        TS_ASSERT_DELTA(testvalue_adams_bashforth,exact_solution,global_error_adams_bashforth);
    }
    
    
    
};

#endif //_TESTABSTRACTIVPODESOLVER_HPP_
