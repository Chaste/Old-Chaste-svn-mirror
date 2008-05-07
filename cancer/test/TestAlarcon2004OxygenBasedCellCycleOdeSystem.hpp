/*

Copyright (C) University of Oxford, 2008

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
#ifndef TESTALARCON2004OXYGENBASEDCELLCYCLEODESYSTEM_HPP_
#define TESTALARCON2004OXYGENBASEDCELLCYCLEODESYSTEM_HPP_

#include <cxxtest/TestSuite.h>

#include <stdio.h>
#include <time.h>
#include <vector>
#include <iostream>

#include "Alarcon2004OxygenBasedCellCycleOdeSystem.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "RungeKuttaFehlbergIvpOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "PetscSetupAndFinalize.hpp"


class TestAlarcon2004OxygenBasedCellCycleOdeSystem : public CxxTest::TestSuite
{
public:

    void TestAlarcon2004Equations()
    {        
        double oxygen_concentration = 1.0;
        Alarcon2004OxygenBasedCellCycleOdeSystem normal_system(oxygen_concentration, HEALTHY);
        Alarcon2004OxygenBasedCellCycleOdeSystem cancer_system(oxygen_concentration, LABELLED);
        
        double time = 0.0;
        std::vector<double> initial_conditions = normal_system.GetInitialConditions();
        
        std::vector<double> normal_derivs(initial_conditions.size());
        std::vector<double> cancer_derivs(initial_conditions.size());
        normal_system.EvaluateYDerivatives(time, initial_conditions, normal_derivs);
        cancer_system.EvaluateYDerivatives(time, initial_conditions, cancer_derivs);
        
        // Test derivatives are correct at t=0 for a normal cell...
        // (figures from MatLab code)
        TS_ASSERT_DELTA(normal_derivs[0], 455.630699088, 1e-5);
        TS_ASSERT_DELTA(normal_derivs[1], 1.83000000000, 1e-5);
        TS_ASSERT_DELTA(normal_derivs[2], 3.00000000000, 1e-5);
        TS_ASSERT_DELTA(normal_derivs[3], 1.50000000000, 1e-5);
        TS_ASSERT_DELTA(normal_derivs[4], -5.4060000000, 1e-5);
        TS_ASSERT_DELTA(normal_derivs[5], 0.00000000000, 1e-5);
        // ...and a cancer cell        
        TS_ASSERT_DELTA(cancer_derivs[0], 455.630699088, 1e-5);
        TS_ASSERT_DELTA(cancer_derivs[1], 1.83600000000, 1e-5);
        TS_ASSERT_DELTA(cancer_derivs[2], 0.42000000000, 1e-5);
        TS_ASSERT_DELTA(cancer_derivs[3], 1.50000000000, 1e-5);
        TS_ASSERT_DELTA(cancer_derivs[4], -5.4060000000, 1e-5);
        TS_ASSERT_DELTA(cancer_derivs[5], 0.00000000000, 1e-5);
        
        // Same thing for a low oxygen concentration 
        // (the usual initial condition for z is zero, so we need to change it 
        // to see any difference)    
        oxygen_concentration = 0.1;
        
        Alarcon2004OxygenBasedCellCycleOdeSystem normal_system2(oxygen_concentration,HEALTHY);
        Alarcon2004OxygenBasedCellCycleOdeSystem cancer_system2(oxygen_concentration,LABELLED);
              
        std::vector<double> normal_derivs2(initial_conditions.size());
        normal_system2.SetInitialConditionsComponent(2, 0.1);
        
        std::vector<double> cancer_derivs2(initial_conditions.size());
        cancer_system2.SetInitialConditionsComponent(2, 0.1);
        
        // initial conditions have changed
        std::vector<double> initial_conditions2 = normal_system2.GetInitialConditions();
                
        normal_system2.EvaluateYDerivatives(time, initial_conditions2, normal_derivs2);
        cancer_system2.EvaluateYDerivatives(time, initial_conditions2, cancer_derivs2);
   
        TS_ASSERT_DELTA(normal_derivs2[0], 455.630699088, 1e-5);
        TS_ASSERT_DELTA(normal_derivs2[1], 1.81500000000, 1e-5);
        TS_ASSERT_DELTA(normal_derivs2[2], 2.94545454545, 1e-5);
        TS_ASSERT_DELTA(normal_derivs2[3], 1.50000000000, 1e-5);
        TS_ASSERT_DELTA(normal_derivs2[4], -5.4060000000, 1e-5);
        TS_ASSERT_DELTA(normal_derivs2[5], 0.00000000000, 1e-5);
        // ...and a cancer cell                
        
        TS_ASSERT_DELTA(cancer_derivs2[0], 455.630699088, 1e-5);
        TS_ASSERT_DELTA(cancer_derivs2[1], 1.82100000000, 1e-5);
        TS_ASSERT_DELTA(cancer_derivs2[2], 0.36545454545, 1e-5);
        TS_ASSERT_DELTA(cancer_derivs2[3], 1.50000000000, 1e-5);
        TS_ASSERT_DELTA(cancer_derivs2[4], -5.4060000000, 1e-5);
        TS_ASSERT_DELTA(cancer_derivs2[5], 0.00000000000, 1e-5);
    }
    
    void TestAlarcon2004Solver() throw(Exception)
    {
        double oxygen_concentration = 1.0;
        Alarcon2004OxygenBasedCellCycleOdeSystem alarcon_system(oxygen_concentration,HEALTHY);
        // Solve system using rk4 solver
        // Matlab's strictest bit uses 0.01 below and relaxes it on flatter bits.
        
        double h_value = 1e-4;
        
        RungeKutta4IvpOdeSolver rk4_solver;
        RungeKuttaFehlbergIvpOdeSolver rkf_solver;
        BackwardEulerIvpOdeSolver back_solver(6);
        
        OdeSolution solutions;
        
        std::vector<double> initial_conditions = alarcon_system.GetInitialConditions();
                
        double start_time, end_time, elapsed_time = 0.0;
        start_time = std::clock();
        solutions = rk4_solver.Solve(&alarcon_system, initial_conditions, 0.0, 10.0, h_value, h_value);
        end_time = std::clock();
        elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
        std::cout <<  "1. Runge-Kutta Elapsed time = " << elapsed_time << "\n";
        
        h_value = 1e-1;
        
        initial_conditions = alarcon_system.GetInitialConditions();
        start_time = std::clock();
        solutions = rkf_solver.Solve(&alarcon_system, initial_conditions, 0.0, 10.0, h_value, 1e-4);
        end_time = std::clock();
        elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
        std::cout <<  "2. Runge-Kutta-Fehlberg Elapsed time = " << elapsed_time << "\n";
        

        // test solutions are OK for a small time increase
        int end = solutions.rGetSolutions().size() - 1;
        
        // tests the simulation is ending at the right time
        TS_ASSERT_DELTA(solutions.rGetTimes()[end] , 9.286356375 , 1e-2);
        
        // proper values found using MatLab ode15s - shocking tolerances though
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0], 0.004000000000000, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1], 0.379221366479055, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2], 0.190488726735972, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][3], 9.962110289977730, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][4], 0.096476600742599, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][5], 1.000000000000000, 1e-3);
    }            
            
};

#endif /*TESTALARCON2004OXYGENBASEDCELLCYCLEODESYSTEM_HPP_*/
