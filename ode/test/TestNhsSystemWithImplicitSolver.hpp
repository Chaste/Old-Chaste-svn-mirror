#ifndef TESTNHSSYSTEMWITHIMPLICITSOLVER_HPP_
#define TESTNHSSYSTEMWITHIMPLICITSOLVER_HPP_

#include <cxxtest/TestSuite.h>
#include "NhsSystemWithImplicitSolver.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "ZeroStimulus.hpp"


class TestNhsSystemWithImplicitSolver : public CxxTest::TestSuite
{
private:
    double GetSampleCaIValue()
    {
        EulerIvpOdeSolver euler_solver;
        ZeroStimulus zero_stimulus;
        LuoRudyIModel1991OdeSystem lr91(&euler_solver, 0.01, &zero_stimulus);
        return lr91.rGetStateVariables()[lr91.GetStateVariableNumberByName("CaI")];
    }        

public:
    void TestSolverSingleTimestep()
    {
        NhsSystemWithImplicitSolver system_with_solver;

        // lam=const, dlamdt not const doesn't make much sense, just for testing purposes
        system_with_solver.SetLambda1AndDerivative(0.5, 0.1);
        double Ca_I = GetSampleCaIValue();
        system_with_solver.SetIntracellularCalciumConcentration(Ca_I);

        // solve system (but don't update state vars yet 
        system_with_solver.SolveDoNotUpdate(0, 0.1, 0.1); // one timestep
        
        NHSCellularMechanicsOdeSystem system_for_euler_solver;
        
        unsigned num_vars = system_with_solver.GetNumberOfStateVariables();
        for(unsigned i=0; i<num_vars; i++)
        {
            // both should be the same (ie initial values)
            TS_ASSERT_DELTA(system_with_solver.rGetStateVariables()[i],  
                            system_for_euler_solver.rGetStateVariables()[i],
                            1e-12);
        }                             

        // solve system with euler
        system_for_euler_solver.SetLambda1AndDerivative(0.5, 0.1); 
        system_for_euler_solver.SetIntracellularCalciumConcentration(Ca_I);
        EulerIvpOdeSolver euler_solver;
        euler_solver.SolveAndUpdateStateVariable(&system_for_euler_solver, 0, 0.1, 0.1);  // one timestep

        // update state vars on implicit system
        system_with_solver.UpdateStateVariables();

        for(unsigned i=0; i<num_vars; i++)
        {
            //std::cout << system_with_solver.rGetStateVariables()[i] << " "
            //          << system_for_euler_solver.rGetStateVariables()[i] << "\n"; 
            
            // we want these within 10% of each other. Note we expect the implicit 
            // solver to be more accurate than the explicit solver, and the timestep 
            // is quite large (as we want non-zero solutions), so can't expect them
            // to be too close.
            TS_ASSERT_DELTA(system_with_solver.rGetStateVariables()[i],  
                            system_for_euler_solver.rGetStateVariables()[i],
                            fabs(system_with_solver.rGetStateVariables()[i]*1e-1));
        }     
    }


    void TestSolverManyTimestepsCompareWithEuler()
    {
        for(unsigned run=0; run<2; run++)
        {
            clock_t ck_start, ck_end;
    
            NhsSystemWithImplicitSolver system_with_solver;
    
            // lam=const, dlamdt not const doesn't make much sense, just for testing purposes
            system_with_solver.SetLambda1AndDerivative(0.5, 0.1);
            double Ca_I = GetSampleCaIValue();
            // bigger Ca_I, so we get some active tension (and so the a few iterations are
            // needed when solving for T_a and z
            system_with_solver.SetIntracellularCalciumConcentration(10*Ca_I);

            // IF THE SECOND RUN, use the implicit-explicit mixed method for z
            if(run==1)
            {
                system_with_solver.UseImplicitExplicitSolveForZ();
            }
    
            // solve system and update
            ck_start = clock();
            system_with_solver.SolveDoNotUpdate(0, 100, 0.01); 
            ck_end = clock();
            double implicit_solve_time = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        
            system_with_solver.UpdateStateVariables();
    
            // solve system with euler
            NHSCellularMechanicsOdeSystem system_for_euler_solver;
            system_for_euler_solver.SetLambda1AndDerivative(0.5, 0.1); 
            system_for_euler_solver.SetIntracellularCalciumConcentration(10*Ca_I);
            EulerIvpOdeSolver euler_solver;

            ck_start = clock();
            euler_solver.SolveAndUpdateStateVariable(&system_for_euler_solver, 0, 100, 0.01); 
            ck_end = clock();
            double explicit_solve_time = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
    
            run==0 ? std::cout<<"\nImplicit vs Explicit\n" : std::cout<<"\nImplicitExplicit vs Explicit\n"; 
            unsigned num_vars = system_with_solver.GetNumberOfStateVariables();
            for(unsigned i=0; i<num_vars; i++)
            {
                std::cout << system_with_solver.rGetStateVariables()[i] << " "
                          << system_for_euler_solver.rGetStateVariables()[i] << "\n"; 
                
                // small timesteps, want these to be very close (1%). and they are
                TS_ASSERT_DELTA(system_with_solver.rGetStateVariables()[i],  
                                system_for_euler_solver.rGetStateVariables()[i],
                                fabs(system_with_solver.rGetStateVariables()[i]*1e-2));
            }     
            std::cout << "TIMES: " << implicit_solve_time << " " << explicit_solve_time << "\n\n";
        }
    }
    
//// test how large a timestep the implicit solver can get away with. needs more study
//    void TestImplicitSolverWithLargeTimeSteps()
//    {
//        NhsSystemWithImplicitSolver system_with_solver;
//        system_with_solver.SetLambda1AndDerivative(0.5, 0.1);
//        system_with_solver.SetIntracellularCalciumConcentration(10*GetSampleCaIValue());
//
//        system_with_solver.SolveDoNotUpdate(0, 100, 0.01); 
//        system_with_solver.UpdateStateVariables();
//
//        NhsSystemWithImplicitSolver system_with_solver2;
//        system_with_solver2.SetLambda1AndDerivative(0.5, 0.1);
//        system_with_solver2.SetIntracellularCalciumConcentration(10*GetSampleCaIValue());
//
//        system_with_solver2.SolveDoNotUpdate(0, 100, 1); 
//        system_with_solver2.UpdateStateVariables();
//
//        unsigned num_vars = system_with_solver.GetNumberOfStateVariables();
//        for(unsigned i=0; i<num_vars; i++)
//        {
//            std::cout << system_with_solver.rGetStateVariables()[i] << " "
//                      << system_with_solver2.rGetStateVariables()[i] << "\n"; 
//            
//            TS_ASSERT_DELTA(system_with_solver.rGetStateVariables()[i],  
//                            system_with_solver2.rGetStateVariables()[i],
//                            fabs(system_with_solver.rGetStateVariables()[i]*1e-2));
//        }
//    }
};

#endif /*TESTNHSSYSTEMWITHIMPLICITSOLVER_HPP_*/
