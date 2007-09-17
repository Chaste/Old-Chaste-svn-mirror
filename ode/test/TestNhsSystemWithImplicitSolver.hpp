#ifndef TESTNHSSYSTEMWITHIMPLICITSOLVER_HPP_
#define TESTNHSSYSTEMWITHIMPLICITSOLVER_HPP_

#include <cxxtest/TestSuite.h>
#include "NhsSystemWithImplicitSolver.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "ZeroStimulus.hpp"


class TestNhsSystemWithImplicitSolver : public CxxTest::TestSuite
{
public:
    void TestSolver()
    {
        NhsSystemWithImplicitSolver system_with_solver;
        system_with_solver.SetLambda1AndDerivative(0.5, 0.1);
        
        unsigned num_vars = system_with_solver.GetNumberOfStateVariables();

        // the following is just to get a realistic Ca_I value
        EulerIvpOdeSolver euler_solver;
        ZeroStimulus zero_stimulus;
        LuoRudyIModel1991OdeSystem lr91(&euler_solver, 0.01, &zero_stimulus);
        unsigned Ca_i_index = lr91.GetStateVariableNumberByName("CaI");
        double Ca_I = lr91.rGetStateVariables()[Ca_i_index];
        
        system_with_solver.SetIntracellularCalciumConcentration(Ca_I);

        system_with_solver.SolveDoNotUpdate(0,0.1,0.1);
        
        NHSCellularMechanicsOdeSystem system_no_solver;
        system_no_solver.SetLambda1AndDerivative(0.5, 0.1);
        system_no_solver.SetIntracellularCalciumConcentration(Ca_I);
        
        for(unsigned i=0; i<num_vars; i++)
        {
            TS_ASSERT_DELTA(system_with_solver.rGetStateVariables()[i],  
                            system_no_solver.rGetStateVariables()[i],
                            1e-12);
        }                             

        euler_solver.SolveAndUpdateStateVariable(&system_no_solver, 0, 0.1, 0.1);
        system_with_solver.UpdateSystem();

        for(unsigned i=0; i<num_vars; i++)
        {
            // we want these within 10% of each other. Note we expect the implicit 
            // solver to be more accurate than the explicit one so can't expect them
            // to be too close.
            TS_ASSERT_DELTA(system_with_solver.rGetStateVariables()[i],  
                            system_no_solver.rGetStateVariables()[i],
                            fabs(system_with_solver.rGetStateVariables()[i]*1e-1));
        }     
    }
};
#endif /*TESTNHSSYSTEMWITHIMPLICITSOLVER_HPP_*/
