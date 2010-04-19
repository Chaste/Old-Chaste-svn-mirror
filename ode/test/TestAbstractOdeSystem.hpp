/*

Copyright (C) University of Oxford, 2005-2010

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


#ifndef _TESTABSTRACTODESYSTEM_HPP_
#define _TESTABSTRACTODESYSTEM_HPP_

// TestAbstractOdeSystem.hpp

#include <cmath>
#include <iostream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "Ode1.hpp"
#include "Ode2.hpp"
#include "Ode3.hpp"
#include "TwoDimOdeSystem.hpp"
#include "VanDerPolOde.hpp"
#include "ParameterisedOde.hpp"
#include "OutputFileHandler.hpp"
// Tolerance for tests
const double tol = 0.01;


class TestAbstractOdeSystem : public CxxTest::TestSuite
{
public:

    void TestOdeSystemOne()
    {
        // Test Ode1 class
        Ode1 ode1;
        // dy
        std::vector<double> dy(1);
        ode1.EvaluateYDerivatives(1.0, ode1.GetInitialConditions(),dy);
        TS_ASSERT_DELTA(dy[0], 1.0, tol);
    }

    void TestOdeSystemTwo()
    {
        Ode2 ode2;
        std::vector<double> dy(1);
        ode2.EvaluateYDerivatives(2.0, ode2.GetInitialConditions(), dy);
        TS_ASSERT_DELTA(dy[0], 8.0, tol);
    }

    void TestOdeSystemThree()
    {
        Ode3 ode3;
        std::vector<double> dy(2);
        ode3.EvaluateYDerivatives(2.0, ode3.GetInitialConditions(), dy);
        TS_ASSERT_DELTA(dy[0], 8.0, tol);
        TS_ASSERT_DELTA(dy[1], 16.0, tol);
    }

    void TestExceptions()
    {
        Ode1 ode;
        TS_ASSERT_EQUALS(ode.GetNumberOfStateVariables(), 1u);
        std::vector<double> v(2);
        v[0] = -1.0;
        v[1] = -2.0;
        TS_ASSERT_THROWS_THIS(ode.SetInitialConditions(v),
                "The number of initial conditions must be that of the number of state variables.");
        TS_ASSERT_THROWS_THIS(ode.SetInitialCondition(2, -3.0),
                "Index is greater than the number of state variables.");
        TS_ASSERT_THROWS_THIS(ode.SetStateVariables(v),
                "The size of the passed in vector must be that of the number of state variables.");
    }

    void TestParameters()
    {
        ParameterisedOde ode;
        boost::shared_ptr<const AbstractOdeSystemInformation> p_info = ode.GetSystemInformation();

        TS_ASSERT_EQUALS(ode.GetParameter(0), 0.0);
        TS_ASSERT_EQUALS(ode.GetNumberOfParameters(), 1u);
        TS_ASSERT_EQUALS(ode.GetNumberOfStateVariables(), 1u);

        double a = 1.0;
        ode.SetParameter(0, a);
        TS_ASSERT_EQUALS(ode.GetParameter(0), a);

        TS_ASSERT_EQUALS(ode.rGetParameterNames()[0], "a");
        TS_ASSERT_EQUALS(ode.rGetParameterUnits()[0], "dimensionless");
        TS_ASSERT_EQUALS(ode.GetParameterIndex("a"), 0u);
        TS_ASSERT_EQUALS(ode.GetParameterUnits(0u), "dimensionless");

        TS_ASSERT_EQUALS(p_info->rGetParameterNames()[0], "a");
        TS_ASSERT_EQUALS(p_info->rGetParameterUnits()[0], "dimensionless");
        TS_ASSERT_EQUALS(p_info->GetParameterIndex("a"), 0u);
        TS_ASSERT_EQUALS(p_info->GetParameterUnits(0u), "dimensionless");

        TS_ASSERT_EQUALS(ode.rGetStateVariableNames()[0], "y");
        TS_ASSERT_EQUALS(ode.rGetStateVariableUnits()[0], "dimensionless");
        TS_ASSERT_EQUALS(ode.GetStateVariableIndex("y"), 0u);
        TS_ASSERT_EQUALS(ode.GetStateVariableUnits(0u), "dimensionless");

        TS_ASSERT_EQUALS(p_info->rGetStateVariableNames()[0], "y");
        TS_ASSERT_EQUALS(p_info->rGetStateVariableUnits()[0], "dimensionless");
        TS_ASSERT_EQUALS(p_info->GetStateVariableIndex("y"), 0u);
        TS_ASSERT_EQUALS(p_info->GetStateVariableUnits(0u), "dimensionless");

        double y = -1.0;
        ode.SetStateVariable(0u, y);
        TS_ASSERT_EQUALS(ode.GetStateVariable(0u), y);

        TS_ASSERT_EQUALS(ode.GetAnyVariableIndex("y"), 0u);
        TS_ASSERT_EQUALS(ode.GetAnyVariableIndex("a"), 1u);
        TS_ASSERT_EQUALS(ode.GetAnyVariable(0u), y);
        TS_ASSERT_EQUALS(ode.GetAnyVariable(1u), a);

        TS_ASSERT_EQUALS(ode.GetAnyVariableUnits(0u), "dimensionless");
        TS_ASSERT_EQUALS(ode.GetAnyVariableUnits(1u), "dimensionless");

        TS_ASSERT_EQUALS(p_info->GetAnyVariableIndex("y"), 0u);
        TS_ASSERT_EQUALS(p_info->GetAnyVariableIndex("a"), 1u);
        TS_ASSERT_EQUALS(p_info->GetAnyVariableUnits(0u), "dimensionless");
        TS_ASSERT_EQUALS(p_info->GetAnyVariableUnits(1u), "dimensionless");

        // Exceptions
        TS_ASSERT_THROWS_THIS(ode.SetParameter(1u, -a), "The index passed in must be less than the number of parameters.");

        TS_ASSERT_THROWS_THIS(ode.GetParameterIndex("b"), "No parameter named 'b'.");
        TS_ASSERT_THROWS_THIS(ode.GetStateVariableIndex("b"), "No state variable named 'b'.");
        TS_ASSERT_THROWS_THIS(ode.GetAnyVariableIndex("b"), "No state variable or parameter named 'b'.");

        TS_ASSERT_THROWS_THIS(ode.GetAnyVariable(2u), "Invalid index passed to GetAnyVariable.");
        TS_ASSERT_THROWS_THIS(ode.GetStateVariable(1u), "The index passed in must be less than the number of state variables.");
        TS_ASSERT_THROWS_THIS(ode.GetParameter(1u), "The index passed in must be less than the number of parameters.");

        TS_ASSERT_THROWS_THIS(ode.GetParameterUnits(1u), "The index passed in must be less than the number of parameters.");
        TS_ASSERT_THROWS_THIS(ode.GetStateVariableUnits(1u), "The index passed in must be less than the number of state variables.");
        TS_ASSERT_THROWS_THIS(ode.GetAnyVariableUnits(2u), "Invalid index passed to GetAnyVariableUnits.");

        TS_ASSERT_THROWS_THIS(p_info->GetParameterIndex("b"), "No parameter named 'b'.");
        TS_ASSERT_THROWS_THIS(p_info->GetStateVariableIndex("b"), "No state variable named 'b'.");
        TS_ASSERT_THROWS_THIS(p_info->GetAnyVariableIndex("b"), "No state variable or parameter named 'b'.");
        TS_ASSERT_THROWS_THIS(p_info->GetParameterUnits(1u), "The index passed in must be less than the number of parameters.");
        TS_ASSERT_THROWS_THIS(p_info->GetStateVariableUnits(1u), "The index passed in must be less than the number of state variables.");
        TS_ASSERT_THROWS_THIS(p_info->GetAnyVariableUnits(2u), "Invalid index passed to GetAnyVariableUnits.");
    }

    void TestSetGetFunctionsInAbstractOdeSystem()
    {
        TwoDimOdeSystem ode;

        std::vector<double> initial_conditions = ode.GetInitialConditions();
        std::vector<double>& r_state_variables = ode.rGetStateVariables();

        TS_ASSERT_DELTA(initial_conditions[0], 1.0, 1e-12);
        TS_ASSERT_DELTA(initial_conditions[1], 2.0, 1e-12);
        TS_ASSERT_DELTA(r_state_variables[0], 3.0, 1e-12);
        TS_ASSERT_DELTA(r_state_variables[1], 4.0, 1e-12);

        std::vector<double> new_initial_conditions;
        new_initial_conditions.push_back(5.0);
        new_initial_conditions.push_back(6.0);

        std::vector<double> new_state_variables;
        new_state_variables.push_back(7.0);
        new_state_variables.push_back(8.0);

        ode.SetInitialConditions(new_initial_conditions);
        ode.SetStateVariables(new_state_variables);

        initial_conditions = ode.GetInitialConditions();

        TS_ASSERT_DELTA(initial_conditions[0], 5.0, 1e-12);
        TS_ASSERT_DELTA(initial_conditions[1], 6.0, 1e-12);
        TS_ASSERT_DELTA(r_state_variables[0], 7.0, 1e-12);
        TS_ASSERT_DELTA(r_state_variables[1], 8.0, 1e-12);

        ode.SetInitialCondition(1, 9.0);
        initial_conditions = ode.GetInitialConditions();
        TS_ASSERT_DELTA(initial_conditions[0], 5.0, 1e-12);
        TS_ASSERT_DELTA(initial_conditions[1], 9.0, 1e-12);

        //Archive the ODE system
        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "ode.arch";
        std::ofstream ofs(archive_filename.c_str());
        boost::archive::text_oarchive output_arch(ofs);

        output_arch <<  static_cast<const TwoDimOdeSystem&>(ode);
        TS_ASSERT_DELTA(r_state_variables[0], 7.0, 1e-12);
        TS_ASSERT_DELTA(r_state_variables[1], 8.0, 1e-12);

        ode.SetStateVariable(0, 2.0);
        ode.SetStateVariable(1, 5.0);

        TS_ASSERT_THROWS_THIS(ode.SetStateVariable(2, 1.0),
                "The index passed in must be less than the number of state variables."); //cover exception

        TS_ASSERT_DELTA(r_state_variables[0], 2.0, 1e-12);
        TS_ASSERT_DELTA(r_state_variables[1], 5.0, 1e-12);
    }
    
    void TestLoadAbstractOdeSystem()
    {
        TwoDimOdeSystem ode;

        TS_ASSERT_EQUALS( ode.GetNumberOfStateVariables(), 2U );

        // Read archive from previous test
        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "ode.arch";

        std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
        boost::archive::text_iarchive input_arch(ifs);

        input_arch >> ode;

        TS_ASSERT_EQUALS( ode.GetNumberOfStateVariables(), 2U );

        std::vector<double>& r_state_variables = ode.rGetStateVariables();

        TS_ASSERT_DELTA(r_state_variables[0], 7.0, 1e-12);
        TS_ASSERT_DELTA(r_state_variables[1], 8.0, 1e-12);
    }

    void TestReadSpecificStateVariable()
    {
        // Create a VanDerPol system
        VanDerPolOde ode_system;

        // get the velocity state variable number
        unsigned var_number = ode_system.GetStateVariableIndex("v");
        TS_ASSERT_EQUALS(var_number, 1u);

        TS_ASSERT_THROWS_THIS(ode_system.GetStateVariableIndex("foo"),
                              "No state variable named 'foo'.");

        TS_ASSERT_EQUALS(ode_system.GetStateVariable(var_number), 10.0);

        TS_ASSERT_EQUALS(ode_system.GetStateVariableUnits(var_number), "m/s");
    }

    // This test is mainly for coverage purposes.
    void TestDumpState()
    {
        // Create a VanDerPol system
        VanDerPolOde ode_system;

        // Dump the state variables
        std::string state = ode_system.DumpState("This is a test.");
        TS_ASSERT_EQUALS(state, "This is a test.\nState:\n\tx:10\n\tv:10\n");

        // Dump user-supplied values
        std::vector<double> rY(2);
        rY[0] = 0.0;
        rY[1] = 1.0;
        state = ode_system.DumpState("Test 2.", rY);
        TS_ASSERT_EQUALS(state, "Test 2.\nState:\n\tx:0\n\tv:1\n");
    }

};

#endif //_TESTABSTRACTODESYSTEM_HPP_
