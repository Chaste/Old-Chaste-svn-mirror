/*

Copyright (C) University of Oxford, 2005-2011

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


#ifndef _TESTABSTRACTCVODESYSTEM_HPP_
#define _TESTABSTRACTCVODESYSTEM_HPP_

#include <cmath>
#include <iostream>
//#include <boost/archive/text_oarchive.hpp>
//#include <boost/archive/text_iarchive.hpp>

#include "Cvode1.hpp"
#include "ParameterisedCvode.hpp"
#include "TwoDimCvodeSystem.hpp"
#include "VectorHelperFunctions.hpp"

double tol = 1e-6;

class TestAbstractCvodeSystem : public CxxTest::TestSuite
{
public:

    void TestOdeSystemOne()
    {
#ifdef CHASTE_CVODE
        // Test Ode1 class
        Cvode1 ode1;
        // dy
        N_Vector y = ode1.GetInitialConditions();
        N_Vector dy = NULL;
        CreateVectorIfEmpty(dy,1);

        ode1.EvaluateYDerivatives(1.0, y, dy);
        TS_ASSERT_DELTA(GetVectorComponent(dy,0), 1.0, tol);
#else
        std::cout << "Cvode is not enabled.\n";
#endif // CHASTE_CVODE
    }

    void TestExceptions()
    {
#ifdef CHASTE_CVODE
        Cvode1 ode;
        TS_ASSERT_EQUALS(ode.GetNumberOfStateVariables(), 1u);

        N_Vector v = NULL;
        CreateVectorIfEmpty(v,2);

        TS_ASSERT_THROWS_THIS(ode.SetDefaultInitialConditions(v),
                "The number of initial conditions must be that of the number of state variables.");
        TS_ASSERT_THROWS_THIS(ode.SetDefaultInitialCondition(2, -3.0),
                "Index is greater than the number of state variables.");
        TS_ASSERT_THROWS_THIS(ode.SetStateVariables(v),
                "The size of the passed in vector must be that of the number of state variables.");
#else
        std::cout << "Cvode is not enabled.\n";
#endif // CHASTE_CVODE
    }

    void TestParameters()
    {
#ifdef CHASTE_CVODE
        ParameterisedCvode ode;
        boost::shared_ptr<const AbstractOdeSystemInformation> p_info = ode.GetSystemInformation();

        TS_ASSERT_EQUALS(ode.GetSystemName(), "ParameterisedCvode");
        TS_ASSERT_EQUALS(p_info->GetSystemName(), "ParameterisedCvode");

        TS_ASSERT_EQUALS(ode.GetParameter(0), 0.0);
        TS_ASSERT_EQUALS(ode.GetNumberOfParameters(), 1u);
        TS_ASSERT_EQUALS(ode.GetNumberOfStateVariables(), 1u);

        double a = 1.0;
        ode.SetParameter(0, a);
        TS_ASSERT_EQUALS(ode.GetParameter(0), a);
        ode.SetParameter("a", a);
        TS_ASSERT_EQUALS(ode.GetParameter(0), a);
        TS_ASSERT_EQUALS(ode.GetParameter("a"), a);

        TS_ASSERT_EQUALS(ode.HasParameter("a"), true);
        TS_ASSERT_EQUALS(ode.rGetParameterNames()[0], "a");
        TS_ASSERT_EQUALS(ode.rGetParameterUnits()[0], "dimensionless");
        TS_ASSERT_EQUALS(ode.GetParameterIndex("a"), 0u);
        TS_ASSERT_EQUALS(ode.GetParameterUnits(0u), "dimensionless");

        TS_ASSERT_EQUALS(p_info->rGetParameterNames()[0], "a");
        TS_ASSERT_EQUALS(p_info->rGetParameterUnits()[0], "dimensionless");
        TS_ASSERT_EQUALS(p_info->GetParameterIndex("a"), 0u);
        TS_ASSERT_EQUALS(p_info->GetParameterUnits(0u), "dimensionless");

        TS_ASSERT_EQUALS(ode.HasStateVariable("y"), true);
        TS_ASSERT_EQUALS(ode.rGetStateVariableNames()[0], "y");
        TS_ASSERT_EQUALS(ode.rGetStateVariableUnits()[0], "dimensionless");
        TS_ASSERT_EQUALS(ode.GetStateVariableIndex("y"), 0u);
        TS_ASSERT_EQUALS(ode.GetStateVariableUnits(0u), "dimensionless");

        TS_ASSERT_EQUALS(p_info->rGetStateVariableNames()[0], "y");
        TS_ASSERT_EQUALS(p_info->rGetStateVariableUnits()[0], "dimensionless");
        TS_ASSERT_EQUALS(p_info->GetStateVariableIndex("y"), 0u);
        TS_ASSERT_EQUALS(p_info->GetStateVariableUnits(0u), "dimensionless");

        double y = -1.0;
        ode.SetAnyVariable(0u, y);
        TS_ASSERT_EQUALS(ode.GetStateVariable("y"), y);
        TS_ASSERT_EQUALS(ode.GetStateVariable(0u), y);

        ode.SetStateVariable("y", y-1);
        TS_ASSERT_EQUALS(ode.GetStateVariable("y"), y-1);
        ode.SetAnyVariable("y", y);
        TS_ASSERT_EQUALS(ode.GetStateVariable("y"), y);


        TS_ASSERT_EQUALS(ode.HasAnyVariable("y"), true);
        TS_ASSERT_EQUALS(ode.HasAnyVariable("a"), true);
        TS_ASSERT_EQUALS(ode.GetAnyVariableIndex("y"), 0u);
        TS_ASSERT_EQUALS(ode.GetAnyVariableIndex("a"), 1u);
        TS_ASSERT_EQUALS(ode.GetAnyVariable(0u), y);
        TS_ASSERT_EQUALS(ode.GetAnyVariable("y"), y);
        TS_ASSERT_EQUALS(ode.GetAnyVariable(1u), a);
        TS_ASSERT_EQUALS(ode.GetAnyVariable("a"), a);
        double new_a = 12345.6789;
        ode.SetAnyVariable(1u, new_a);
        TS_ASSERT_EQUALS(ode.GetAnyVariable(1u), new_a);

        TS_ASSERT_EQUALS(ode.GetAnyVariableUnits(0u), "dimensionless");
        TS_ASSERT_EQUALS(ode.GetAnyVariableUnits(1u), "dimensionless");
        TS_ASSERT_EQUALS(ode.GetAnyVariableUnits("y"), "dimensionless");
        TS_ASSERT_EQUALS(ode.GetAnyVariableUnits("a"), "dimensionless");

        TS_ASSERT_EQUALS(p_info->GetAnyVariableIndex("y"), 0u);
        TS_ASSERT_EQUALS(p_info->GetAnyVariableIndex("a"), 1u);
        TS_ASSERT_EQUALS(p_info->GetAnyVariableUnits(0u), "dimensionless");
        TS_ASSERT_EQUALS(p_info->GetAnyVariableUnits(1u), "dimensionless");

        // Exceptions
        TS_ASSERT_THROWS_THIS(ode.SetParameter(1u, -a), "The index passed in must be less than the number of parameters.");

        TS_ASSERT_EQUALS(ode.HasParameter("b"), false);
        TS_ASSERT_THROWS_THIS(ode.GetParameterIndex("b"), "No parameter named 'b'.");
        TS_ASSERT_THROWS_THIS(ode.GetParameter("b"), "No parameter named 'b'.");
        TS_ASSERT_EQUALS(ode.HasStateVariable("b"), false);
        TS_ASSERT_THROWS_THIS(ode.GetStateVariableIndex("b"), "No state variable named 'b'.");
        TS_ASSERT_EQUALS(ode.HasAnyVariable("b"), false);
        TS_ASSERT_THROWS_THIS(ode.GetAnyVariableIndex("b"), "No state variable, parameter, or derived quantity named 'b'.");
        TS_ASSERT_THROWS_THIS(ode.SetAnyVariable(2u, 0.0), "Cannot set the value of a derived quantity, or invalid index.");

        TS_ASSERT_THROWS_THIS(ode.GetAnyVariable(3u), "Invalid index passed to GetAnyVariable.");
        TS_ASSERT_THROWS_THIS(ode.GetStateVariable(1u), "The index passed in must be less than the number of state variables.");
        TS_ASSERT_THROWS_THIS(ode.GetParameter(1u), "The index passed in must be less than the number of parameters.");

        TS_ASSERT_THROWS_THIS(ode.GetParameterUnits(1u), "The index passed in must be less than the number of parameters.");
        TS_ASSERT_THROWS_THIS(ode.GetStateVariableUnits(1u), "The index passed in must be less than the number of state variables.");
        TS_ASSERT_THROWS_THIS(ode.GetAnyVariableUnits(3u), "Invalid index passed to GetAnyVariableUnits.");

        TS_ASSERT_EQUALS(p_info->HasParameter("b"), false);
        TS_ASSERT_THROWS_THIS(p_info->GetParameterIndex("b"), "No parameter named 'b'.");
        TS_ASSERT_EQUALS(p_info->HasStateVariable("b"), false);
        TS_ASSERT_THROWS_THIS(p_info->GetStateVariableIndex("b"), "No state variable named 'b'.");
        TS_ASSERT_EQUALS(p_info->HasAnyVariable("b"), false);
        TS_ASSERT_THROWS_THIS(p_info->GetAnyVariableIndex("b"), "No state variable, parameter, or derived quantity named 'b'.");
        TS_ASSERT_THROWS_THIS(p_info->GetParameterUnits(1u), "The index passed in must be less than the number of parameters.");
        TS_ASSERT_THROWS_THIS(p_info->GetStateVariableUnits(1u), "The index passed in must be less than the number of state variables.");
        TS_ASSERT_THROWS_THIS(p_info->GetAnyVariableUnits(3u), "Invalid index passed to GetAnyVariableUnits.");
#else
        std::cout << "Cvode is not enabled.\n";
#endif // CHASTE_CVODE
    }

    void TestAttributes() throw (Exception)
    {
#ifdef CHASTE_CVODE
        ParameterisedCvode ode;
        TS_ASSERT_EQUALS(ode.GetNumberOfAttributes(), 1u);
        TS_ASSERT(ode.HasAttribute("attr"));
        TS_ASSERT_DELTA(ode.GetAttribute("attr"), 1.1, 1e-12);
        TS_ASSERT(!ode.HasAttribute("missing"));
        TS_ASSERT_THROWS_THIS(ode.GetAttribute("missing"), "No attribute 'missing' found.");

        boost::shared_ptr<const AbstractOdeSystemInformation> p_info = ode.GetSystemInformation();
        TS_ASSERT_EQUALS(p_info->GetNumberOfAttributes(), 1u);
        TS_ASSERT(p_info->HasAttribute("attr"));
        TS_ASSERT_DELTA(p_info->GetAttribute("attr"), 1.1, 1e-12);
        TS_ASSERT(!p_info->HasAttribute("missing"));
        TS_ASSERT_THROWS_THIS(p_info->GetAttribute("missing"), "No attribute 'missing' found.");
#else
        std::cout << "Cvode is not enabled.\n";
#endif // CHASTE_CVODE
    }

    void TestDerivedQuantities() throw (Exception)
    {
#ifdef CHASTE_CVODE
        ParameterisedCvode ode;
        boost::shared_ptr<const AbstractOdeSystemInformation> p_info = ode.GetSystemInformation();

        TS_ASSERT_EQUALS(ode.GetNumberOfDerivedQuantities(), 1u);
        TS_ASSERT_EQUALS(ode.HasDerivedQuantity("2a_plus_y"), true);
        TS_ASSERT_EQUALS(ode.HasDerivedQuantity("Not_there"), false);
        TS_ASSERT_EQUALS(ode.GetDerivedQuantityIndex("2a_plus_y"), 0u);
        TS_ASSERT_EQUALS(ode.GetDerivedQuantityUnits(0u), "dimensionless");
        TS_ASSERT_EQUALS(ode.rGetDerivedQuantityNames().size(), 1u);
        TS_ASSERT_EQUALS(ode.rGetDerivedQuantityUnits().size(), 1u);
        TS_ASSERT_EQUALS(ode.rGetDerivedQuantityNames()[0], "2a_plus_y");
        TS_ASSERT_EQUALS(ode.rGetDerivedQuantityUnits()[0], "dimensionless");

        TS_ASSERT_EQUALS(p_info->HasDerivedQuantity("2a_plus_y"), true);
        TS_ASSERT_EQUALS(p_info->HasDerivedQuantity("Not_there"), false);
        TS_ASSERT_EQUALS(p_info->GetDerivedQuantityIndex("2a_plus_y"), 0u);
        TS_ASSERT_EQUALS(p_info->GetDerivedQuantityUnits(0u), "dimensionless");
        TS_ASSERT_EQUALS(p_info->rGetDerivedQuantityNames().size(), 1u);
        TS_ASSERT_EQUALS(p_info->rGetDerivedQuantityUnits().size(), 1u);
        TS_ASSERT_EQUALS(p_info->rGetDerivedQuantityNames()[0], "2a_plus_y");
        TS_ASSERT_EQUALS(p_info->rGetDerivedQuantityUnits()[0], "dimensionless");

        TS_ASSERT_EQUALS(ode.HasAnyVariable("2a_plus_y"), true);
        TS_ASSERT_EQUALS(ode.HasAnyVariable("Not_there"), false);
        TS_ASSERT_EQUALS(ode.GetAnyVariableIndex("2a_plus_y"), 2u);
        TS_ASSERT_EQUALS(ode.GetAnyVariableUnits(2u), "dimensionless");
        TS_ASSERT_EQUALS(p_info->GetAnyVariableIndex("2a_plus_y"), 2u);
        TS_ASSERT_EQUALS(p_info->GetAnyVariableUnits(2u), "dimensionless");

        N_Vector derived = ode.ComputeDerivedQuantitiesFromCurrentState(0.0);
        double a = ode.GetParameter(0);
        TS_ASSERT_EQUALS(a, 0.0);
        TS_ASSERT_DELTA(GetVectorComponent(derived,0), 2*a, 1e-4);
        TS_ASSERT_DELTA(ode.GetAnyVariable(2u, 0.0), 2*a, 1e-4);
        TS_ASSERT_DELTA(ode.GetAnyVariable(2u, 0.0, &derived), 2*a, 1e-4);
        a = 1.0;
        ode.SetParameter(0, a);
        derived = ode.ComputeDerivedQuantities(0.0, ode.GetInitialConditions());
        TS_ASSERT_DELTA(GetVectorComponent(derived,0), 2*a, 1e-4);
        double y = 10.0;
        ode.SetStateVariable(0, y);
        derived = ode.ComputeDerivedQuantitiesFromCurrentState(0.0);
        TS_ASSERT_DELTA(GetVectorComponent(derived,0), 2*a+y, 1e-4);
        TS_ASSERT_DELTA(ode.GetAnyVariable(2u, 1.0/* ignored for this ODE */), 2*a+y, 1e-4);

        // Exceptions
        TS_ASSERT_THROWS_THIS(ode.GetDerivedQuantityIndex("Missing"), "No derived quantity named 'Missing'.");
        TS_ASSERT_THROWS_THIS(ode.GetDerivedQuantityUnits(1u), "The index passed in must be less than the number of derived quantities.");

        TwoDimCvodeSystem ode2;
        TS_ASSERT_THROWS_THIS(ode2.ComputeDerivedQuantitiesFromCurrentState(0.0),
                              "This ODE system does not define derived quantities.");
        N_Vector doesnt_matter;
        TS_ASSERT_THROWS_THIS(ode2.ComputeDerivedQuantities(0.0, doesnt_matter),
                              "This ODE system does not define derived quantities.");
#else
        std::cout << "Cvode is not enabled.\n";
#endif // CHASTE_CVODE
    }

    // This test is mainly for coverage purposes.
    void TestDumpState()
    {
#ifdef CHASTE_CVODE
        // Create a two variable system
        TwoDimCvodeSystem ode_system;

        // Dump the state variables
        std::string state = ode_system.DumpState("This is a test.");
        TS_ASSERT_EQUALS(state, "This is a test.\nState:\n\tVariable_1:1\n\tVariable_2:2\n");

        // Dump user-supplied values
        N_Vector rY = NULL;
        CreateVectorIfEmpty(rY,2);
        SetVectorComponent(rY, 0, 0.0);
        SetVectorComponent(rY, 1, 1.0);

        state = ode_system.DumpState("Test 2.", rY);
        TS_ASSERT_EQUALS(state, "Test 2.\nState:\n\tVariable_1:0\n\tVariable_2:1\n");
#else
        std::cout << "Cvode is not enabled.\n";
#endif // CHASTE_CVODE
    }
};

#endif //_TESTABSTRACTODESYSTEM_HPP_
