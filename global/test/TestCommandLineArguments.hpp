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


#ifndef TESTCOMMANDLINEARGUMENTS_HPP_
#define TESTCOMMANDLINEARGUMENTS_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>
#include "CommandLineArguments.hpp"

class TestCommandLineArguments : public CxxTest::TestSuite
{
public:
    void TestCommandLineArgumentsSingleton() throw(Exception)
    {
        // test that argc and argv are populated
        int argc = *(CommandLineArguments::Instance()->p_argc);
        TS_ASSERT_LESS_THAN(0, argc); // argc should always be 1 or greater
        
        // argv[0] will be equal to global/build/debug/TestCommandLineArgumentsRunner
        // or global/build/optimised/TestCommandLineArgumentsRunner, etc
        char** argv = *(CommandLineArguments::Instance()->p_argv);
        assert(argv != NULL);
        std::string arg_as_string(argv[0]);
        std::string final_part_of_string = arg_as_string.substr(arg_as_string.length()-30,arg_as_string.length());
        TS_ASSERT_EQUALS("TestCommandLineArgumentsRunner",final_part_of_string);
        
        // Now test OptionExists() and GetValueCorrespondingToOption()
        //
        // The following tests would require the following arguments to be passed
        // in: 
        // ./global/build/debug/TestCommandLineArgumentsRunner -myoption -myintval 24 -mydoubleval 3.14
        // 
        // To test the methods we overwrite the arg_c and arg_v contained in the 
        // singleton with the arguments that were needed.        
        int new_argc = 6;
        char new_argv0[] = "..";
        char new_argv1[] = "-myoption";
        char new_argv2[] = "-myintval";
        char new_argv3[] = "24";
        char new_argv4[] = "-mydoubleval";
        char new_argv5[] = "3.14";
        
        char** new_argv = new char*[6];
        new_argv[0] = new_argv0;
        new_argv[1] = new_argv1;
        new_argv[2] = new_argv2;
        new_argv[3] = new_argv3;
        new_argv[4] = new_argv4;
        new_argv[5] = new_argv5;

        // (save the real args to be restored at the end)
        int* p_real_argc = CommandLineArguments::Instance()->p_argc;
        char*** p_real_argv = CommandLineArguments::Instance()->p_argv;

        // overwrite the args
        CommandLineArguments::Instance()->p_argc = &new_argc;
        CommandLineArguments::Instance()->p_argv = &new_argv;
        
        // test OptionExists()
        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-myoption"));        
        TS_ASSERT( ! CommandLineArguments::Instance()->OptionExists("-asddsgijdfgokgfgurgher"));
        
        // test GetValueCorrespondingToOption()
        char* val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("-myintval");
        unsigned i = atol(val);
        TS_ASSERT_EQUALS(i, 24u);

        val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("-mydoubleval");
        double x = atof(val);
        TS_ASSERT_EQUALS(x, 3.14);

        // test exceptions in GetValueCorrespondingToOption()
        TS_ASSERT_THROWS_CONTAINS(CommandLineArguments::Instance()->GetValueCorrespondingToOption("-rwesdb"), "does not exist");

        new_argc = 5;
        TS_ASSERT_THROWS_CONTAINS(CommandLineArguments::Instance()->GetValueCorrespondingToOption("-mydoubleval"), "No value given after");

        delete new_argv;

        // (restore the real args)
        CommandLineArguments::Instance()->p_argc = p_real_argc;
        CommandLineArguments::Instance()->p_argv = p_real_argv;
    }
};



#endif /*TESTCOMMANDLINEARGUMENTS_HPP_*/
