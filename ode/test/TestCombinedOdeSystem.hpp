/*

Copyright (C) University of Oxford, 2005-2009

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
#ifndef TESTCOMBINEDODESYSTEM_HPP_
#define TESTCOMBINEDODESYSTEM_HPP_

#include <cxxtest/TestSuite.h>

#include "CombinedOdeSystem.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSystemInformation.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "OdeSolution.hpp"

/**
 * dy/dt = x, y(0) = 0. Here x is a parameter.
 */
class SimpleOde1 : public AbstractOdeSystem
{
public:
    SimpleOde1() : AbstractOdeSystem(1) // 1 here is the number of variables
    {
        mpSystemInfo = OdeSystemInformation<SimpleOde1>::Instance();
        SetStateVariables(GetInitialConditions());
        mParameters.resize(1);
    }

    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
    {
        rDY[0] = mParameters[0];
    }
};

template<>
void OdeSystemInformation<SimpleOde1>::Initialise()
{
    this->mVariableNames.push_back("Variable y");
    this->mVariableUnits.push_back("Units y");
    this->mInitialConditions.push_back(0.0);
    
    this->mParameterNames.push_back("Variable x");
    this->mParameterUnits.push_back("Units x");
    
    this->mInitialised = true;
}


/**
 * dx/dt = -y, x(0) = 1. Here y is a parameter.
 */
class SimpleOde2 : public AbstractOdeSystem
{
public:
    SimpleOde2() : AbstractOdeSystem(1) // 1 here is the number of variables
    {
        mpSystemInfo = OdeSystemInformation<SimpleOde2>::Instance();
        SetStateVariables(GetInitialConditions());
        mParameters.resize(1);
    }

    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
    {
        rDY[0] = -mParameters[0];
    }
};

template<>
void OdeSystemInformation<SimpleOde2>::Initialise()
{
    this->mVariableNames.push_back("Variable x");
    this->mVariableUnits.push_back("Units x");
    this->mInitialConditions.push_back(1.0);
    
    this->mParameterNames.push_back("Variable y");
    this->mParameterUnits.push_back("Units y");
    
    this->mInitialised = true;
}



/**
 * The following classes are used in the solution of
 * x'=x-y+z, y'=y-z, z'=2y-z
 * starting at (x,y,z)=(0,1,0)
 * Analytic solution is x=-sin(t), y=sin(t)+cos(t), z=2sin(t)
 */

/**
 * dx/dt = x-y+z, x(0) = 0. Here y & z are parameters.
 */
class SimpleOde3 : public AbstractOdeSystem
{
public:
    SimpleOde3() : AbstractOdeSystem(1) // 1 here is the number of variables
    {
        mpSystemInfo = OdeSystemInformation<SimpleOde3>::Instance();
        SetStateVariables(GetInitialConditions());
        mParameters.resize(2);
    }

    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
    {
        rDY[0] = rY[0] - mParameters[0] + mParameters[1];
    }
};

template<>
void OdeSystemInformation<SimpleOde3>::Initialise()
{
    this->mVariableNames.push_back("x");
    this->mVariableUnits.push_back("dimensionless_x");
    this->mInitialConditions.push_back(0.0);
    
    this->mInitialised = true;
}

/**
 * dy/dt = y-z, y(0) = 1. Here z is a parameter.
 */
class SimpleOde4 : public AbstractOdeSystem
{
public:
    SimpleOde4() : AbstractOdeSystem(1) // 1 here is the number of variables
    {
        mpSystemInfo = OdeSystemInformation<SimpleOde4>::Instance();
        SetStateVariables(GetInitialConditions());
        mParameters.resize(1);
    }

    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
    {
        rDY[0] = rY[0] - mParameters[0];
    }
};

template<>
void OdeSystemInformation<SimpleOde4>::Initialise()
{
    this->mVariableNames.push_back("y");
    this->mVariableUnits.push_back("dimensionless_y");
    this->mInitialConditions.push_back(1.0);
    
    this->mInitialised = true;
}


/**
 * dz/dt = 2y-z, z(0) = 0. Here y is a parameter.
 */
class SimpleOde5 : public AbstractOdeSystem
{
public:
    SimpleOde5() : AbstractOdeSystem(1) // 1 here is the number of variables
    {
        mpSystemInfo = OdeSystemInformation<SimpleOde5>::Instance();
        SetStateVariables(GetInitialConditions());
        mParameters.resize(1);
    }

    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
    {
        rDY[0] = 2.0*mParameters[0] - rY[0];
    }
};

template<>
void OdeSystemInformation<SimpleOde5>::Initialise()
{
    this->mVariableNames.push_back("z");
    this->mVariableUnits.push_back("dimensionless_z");
    this->mInitialConditions.push_back(0.0);
    
    this->mInitialised = true;
}


/**
 * Double system (combination of SimpleOde4, for y' and SimpleOde5, for z') to be used with SimpleOde3
 * dy/dt = y-z, y(0) = 1.
 * dz/dt = 2y-z, z(0) = 0.
 */
class SimpleOde6 : public AbstractOdeSystem
{
public:
    SimpleOde6() : AbstractOdeSystem(2) // 2 here is the number of variables
    {
        mpSystemInfo = OdeSystemInformation<SimpleOde6>::Instance();
        SetStateVariables(GetInitialConditions());
    }

    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
    {
        rDY[0] = rY[0] - rY[1];        //y'=y-z
        rDY[1] = 2.0*rY[0] - rY[1];    //z'=2y-z
    }
};

template<>
void OdeSystemInformation<SimpleOde6>::Initialise()
{
    this->mVariableNames.push_back("y");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(1.0);

    this->mVariableNames.push_back("z");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.0);
    
    this->mInitialised = true;
}


class TestCombinedOdeSystem : public CxxTest::TestSuite
{
public:
    
    void TestSimpleCombinedOdeSystem() throw (Exception)
    {
        // Create two ODE systems
        SimpleOde1 ode_for_y; // dy/dt = x
        SimpleOde2 ode_for_x; // dx/dt = -y

        std::vector<AbstractOdeSystem*> ode_systems;
        ode_systems.push_back(&ode_for_y);
        ode_systems.push_back(&ode_for_x);

        // Create combined ODE system
        CombinedOdeSystem combined_ode_system(ode_systems);

        // Tell the combined ODE system which state variables in the first ODE system
        // correspond to which parameters in the second ODE system...
        std::map<unsigned, unsigned> variable_parameter_map;
        variable_parameter_map[0] = 0;

        combined_ode_system.Configure(variable_parameter_map, &ode_for_y, &ode_for_x);

        // ...and vice versa (we can re-use the map in this case)
        combined_ode_system.Configure(variable_parameter_map, &ode_for_x, &ode_for_y);

        // Test number of state variables
        unsigned num_variables = combined_ode_system.GetNumberOfStateVariables();
        TS_ASSERT_EQUALS(num_variables, 2u);
        
        // Combined system has no parameters
        TS_ASSERT_EQUALS(combined_ode_system.GetNumberOfParameters(), 0u);
        TS_ASSERT_EQUALS(combined_ode_system.rGetParameterNames().size(), 0u);
        
        // Test initial conditions
        std::vector<double> initial_conditions = combined_ode_system.GetInitialConditions();
        TS_ASSERT_DELTA(initial_conditions[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(initial_conditions[1], 1.0, 1e-12);
        // Test variable names & units
        const std::vector<std::string>& r_names = combined_ode_system.rGetVariableNames();
        TS_ASSERT_EQUALS(r_names[0], ode_for_y.rGetVariableNames()[0]);
        TS_ASSERT_EQUALS(r_names[1], ode_for_x.rGetVariableNames()[0]);
        const std::vector<std::string>& r_units = combined_ode_system.rGetVariableUnits();
        TS_ASSERT_EQUALS(r_units[0], ode_for_y.rGetVariableUnits()[0]);
        TS_ASSERT_EQUALS(r_units[1], ode_for_x.rGetVariableUnits()[0]);
        
        // Test solving the combined system.
        // This is dy/dt = x, dx/dt = -y, y(0) = 0, x(0) = 1.
        // The analytic solution is y = sin(t), x = cos(t).
        EulerIvpOdeSolver solver;
        OdeSolution solutions;
        double h = 0.01;
        std::vector<double> inits = combined_ode_system.GetInitialConditions(); 
        solutions = solver.Solve(&combined_ode_system, inits, 0.0, 2.0, h, h);
        double global_error = 0.5*(exp(2)-1)*h;
        TS_ASSERT_DELTA(solutions.rGetSolutions().back()[0], sin(2), global_error);
        TS_ASSERT_DELTA(solutions.rGetSolutions().back()[1], cos(2), global_error);
        
        // Check that if we create the same combination, we get the same information object
        boost::shared_ptr<const AbstractOdeSystemInformation> info1 = combined_ode_system.GetSystemInformation();
        CombinedOdeSystem combined_ode_system2(ode_systems);
        boost::shared_ptr<const AbstractOdeSystemInformation> info2 = combined_ode_system2.GetSystemInformation();
        TS_ASSERT_EQUALS(info1, info2);
    }
    
    void TestSimpleSystemWithOrderSwapped()
    {
        // The solution should be the same, but we'll have to construct a new CombinedOdeSystemInformation
        // object, because the order of subsystems has changed.
        
        SimpleOde1 ode_for_y; // dy/dt = x
        SimpleOde2 ode_for_x; // dx/dt = -y

        std::vector<AbstractOdeSystem*> ode_systems;
        ode_systems.push_back(&ode_for_x);
        ode_systems.push_back(&ode_for_y);

        // Create combined ODE system
        CombinedOdeSystem combined_ode_system(ode_systems);

        // Tell the combined ODE system which state variables in the first ODE system
        // correspond to which parameters in the second ODE system...
        std::map<unsigned, unsigned> variable_parameter_map;
        variable_parameter_map[0] = 0;

        combined_ode_system.Configure(variable_parameter_map, &ode_for_y, &ode_for_x);

        // ...and vice versa (we can re-use the map in this case)
        combined_ode_system.Configure(variable_parameter_map, &ode_for_x, &ode_for_y);

        // Test solving the combined system.
        // This is dy/dt = x, dx/dt = -y, y(0) = 0, x(0) = 1.
        // The analytic solution is y = sin(t), x = cos(t).
        EulerIvpOdeSolver solver;
        OdeSolution solutions;
        double h = 0.01;
        std::vector<double> inits = combined_ode_system.GetInitialConditions(); 
        solutions = solver.Solve(&combined_ode_system, inits, 0.0, 2.0, h, h);
        double global_error = 0.5*(exp(2)-1)*h;
        TS_ASSERT_DELTA(solutions.rGetSolutions().back()[1], sin(2), global_error);
        TS_ASSERT_DELTA(solutions.rGetSolutions().back()[0], cos(2), global_error);
    }
    
    void TestWithThreeVariables()
    {
        SimpleOde3 ode_for_x; // dx/dt = x -y +z
        SimpleOde4 ode_for_y; // dy/dt = y-z
        SimpleOde5 ode_for_z; // dz/dt = 2y-z
        
        std::vector<AbstractOdeSystem*> ode_systems;
        ode_systems.push_back(&ode_for_x);
        ode_systems.push_back(&ode_for_y);
        ode_systems.push_back(&ode_for_z);

        // Create combined ODE system
        CombinedOdeSystem combined_ode_system(ode_systems);

        // Tell the combined ODE system which state variables in the first ODE system
        // correspond to which parameters in the second ODE system...
        std::map<unsigned, unsigned> variable_parameter_map_yx;
        variable_parameter_map_yx[0] = 0; //y in the y-ODE appears in the x-ODE
        combined_ode_system.Configure(variable_parameter_map_yx, &ode_for_y, &ode_for_x);

        std::map<unsigned, unsigned> variable_parameter_map_zx;
        variable_parameter_map_zx[0] = 1; //z in the z-ODE appears in the x-ODE
        combined_ode_system.Configure(variable_parameter_map_zx, &ode_for_z, &ode_for_x);
        
        //Reuse the variable_parameter_map_yx
        //y in the y-ODE appears in the z-ODE
        combined_ode_system.Configure(variable_parameter_map_yx, &ode_for_y, &ode_for_z);
        //z in the z-ODE appears in the y-ODE
        combined_ode_system.Configure(variable_parameter_map_yx, &ode_for_z, &ode_for_y);
        

        // Test number of state variables
        unsigned num_variables = combined_ode_system.GetNumberOfStateVariables();
        TS_ASSERT_EQUALS(num_variables, 3u);
        
        // Combined system has no parameters
        TS_ASSERT_EQUALS(combined_ode_system.GetNumberOfParameters(), 0u);
        TS_ASSERT_EQUALS(combined_ode_system.rGetParameterNames().size(), 0u);
        
        // Test initial conditions
        std::vector<double> initial_conditions = combined_ode_system.GetInitialConditions();
        TS_ASSERT_DELTA(initial_conditions[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(initial_conditions[1], 1.0, 1e-12);
        TS_ASSERT_DELTA(initial_conditions[2], 0.0, 1e-12);
        // Test variable names & units
        const std::vector<std::string>& r_names = combined_ode_system.rGetVariableNames();
        TS_ASSERT_EQUALS(r_names[0], ode_for_x.rGetVariableNames()[0]);
        TS_ASSERT_EQUALS(r_names[1], ode_for_y.rGetVariableNames()[0]);
        TS_ASSERT_EQUALS(r_names[2], ode_for_z.rGetVariableNames()[0]);
        const std::vector<std::string>& r_units = combined_ode_system.rGetVariableUnits();
        TS_ASSERT_EQUALS(r_units[0], ode_for_x.rGetVariableUnits()[0]);
        TS_ASSERT_EQUALS(r_units[1], ode_for_y.rGetVariableUnits()[0]);
        TS_ASSERT_EQUALS(r_units[2], ode_for_z.rGetVariableUnits()[0]);
        TS_ASSERT_EQUALS(r_units[0], "dimensionless_x");
        TS_ASSERT_EQUALS(r_units[1], "dimensionless_y");
        TS_ASSERT_EQUALS(r_units[2], "dimensionless_z");
        
        // x'=x-y+z, y'=y-z, z'=2y-z
        // starting at (x,y,z)=(0,1,0)
        // Analytic solution is x=-sin(t), y=sin(t)+cos(t), z=2sin(t)
        EulerIvpOdeSolver solver;
        OdeSolution solutions;
        double h = 0.01;
        std::vector<double> inits = combined_ode_system.GetInitialConditions(); 
        solutions = solver.Solve(&combined_ode_system, inits, 0.0, 2.0, h, h);
        double global_error = 0.5*(exp(2)-1)*h;
        TS_ASSERT_DELTA(solutions.rGetSolutions().back()[0], -sin(2), global_error);
        TS_ASSERT_DELTA(solutions.rGetSolutions().back()[1], sin(2)+cos(2), global_error);
        TS_ASSERT_DELTA(solutions.rGetSolutions().back()[2], 2.0*sin(2), global_error);
    }

    void TestWithThreeVariablesTwoSystems()
    {
        SimpleOde3 ode_for_x; // dx/dt = x -y +z
        SimpleOde6 ode_for_yz; // dy/dt = y-z  and dz/dt = 2y-z
        
        std::vector<AbstractOdeSystem*> ode_systems;
        ode_systems.push_back(&ode_for_x);
        ode_systems.push_back(&ode_for_yz);

        // Create combined ODE system
        CombinedOdeSystem combined_ode_system(ode_systems);

        // Tell the combined ODE system which state variables in the first ODE system
        // correspond to which parameters in the second ODE system...
        std::map<unsigned, unsigned> variable_parameter_map;
        variable_parameter_map[0] = 0; //y in the yz-ODE appears in the x-ODE
        variable_parameter_map[1] = 1; //z in the yz-ODE appears in the x-ODE
        combined_ode_system.Configure(variable_parameter_map, &ode_for_yz, &ode_for_x);

        // Test number of state variables
        unsigned num_variables = combined_ode_system.GetNumberOfStateVariables();
        TS_ASSERT_EQUALS(num_variables, 3u);

        // Combined system has no parameters
        TS_ASSERT_EQUALS(combined_ode_system.GetNumberOfParameters(), 0u);
        TS_ASSERT_EQUALS(combined_ode_system.rGetParameterNames().size(), 0u);
        
        // Test initial conditions
        std::vector<double> initial_conditions = combined_ode_system.GetInitialConditions();
        TS_ASSERT_DELTA(initial_conditions[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(initial_conditions[1], 1.0, 1e-12);
        TS_ASSERT_DELTA(initial_conditions[2], 0.0, 1e-12);
        // Test variable names & units
        const std::vector<std::string>& r_names = combined_ode_system.rGetVariableNames();
        TS_ASSERT_EQUALS(r_names[0], ode_for_x.rGetVariableNames()[0]);
        TS_ASSERT_EQUALS(r_names[1], ode_for_yz.rGetVariableNames()[0]);
        TS_ASSERT_EQUALS(r_names[2], ode_for_yz.rGetVariableNames()[1]);
        
        // x'=x-y+z, y'=y-z, z'=2y-z
        // starting at (x,y,z)=(0,1,0)
        // Analytic solution is x=-sin(t), y=sin(t)+cos(t), z=2sin(t)
        EulerIvpOdeSolver solver;
        OdeSolution solutions;
        double h = 0.01;
        std::vector<double> inits = combined_ode_system.GetInitialConditions(); 
        solutions = solver.Solve(&combined_ode_system, inits, 0.0, 2.0, h, h);
        double global_error = 0.5*(exp(2)-1)*h;
        TS_ASSERT_DELTA(solutions.rGetSolutions().back()[0], -sin(2), global_error);
        TS_ASSERT_DELTA(solutions.rGetSolutions().back()[1], sin(2)+cos(2), global_error);
        TS_ASSERT_DELTA(solutions.rGetSolutions().back()[2], 2.0*sin(2), global_error);
    }
};

#endif /*TESTCOMBINEDODESYSTEM_HPP_*/
