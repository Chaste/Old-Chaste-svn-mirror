#ifndef TESTIMPLICITNHSCELLULARMECHSYSTEMSOLVER_HPP_
#define TESTIMPLICITNHSCELLULARMECHSYSTEMSOLVER_HPP_

#include <cxxtest/TestSuite.h>
#include "ImplicitNhsCellularMechSystemSolver.hpp"
#include "EulerIvpOdeSolver.hpp"

class TestImplicitNhsCellularMechSystemSolver : public CxxTest::TestSuite
{
public:
    void TestSolver()
    {
        NHSCellularMechanicsOdeSystem system;
        ImplicitNhsCellularMechSystemSolver solver(system);

//set CaI

        solver.SolveDoNotUpdate(0,1,0.01);
        
        NHSCellularMechanicsOdeSystem system2;
        for(unsigned i=0; i<system.GetNumberOfStateVariables(); i++)
        {
            std::cout << system.rGetStateVariables()[i] << "\n";
            TS_ASSERT_DELTA(system.rGetStateVariables()[i],  
                            system2.rGetStateVariables()[i],
                            1e-12);
        }                             

        EulerIvpOdeSolver euler_solver;
        euler_solver.SolveAndUpdateStateVariable(&system2, 0, 1, 0.01);
        solver.UpdateSystem();

        for(unsigned i=0; i<system.GetNumberOfStateVariables(); i++)
        {
            TS_ASSERT_DELTA(system.rGetStateVariables()[i],  
                            system2.rGetStateVariables()[i],
                            system.rGetStateVariables()[i]*1e-3);
        }     
    }
};
#endif /*TESTIMPLICITNHSCELLULARMECHSYSTEMSOLVER_HPP_*/
