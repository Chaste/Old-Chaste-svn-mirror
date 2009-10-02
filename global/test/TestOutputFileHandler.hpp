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


#ifndef TESTOUTPUTFILEHANDLER_HPP_
#define TESTOUTPUTFILEHANDLER_HPP_

#include <cxxtest/TestSuite.h>
#include <string>
#include <fstream>
#include <sstream>
#include <unistd.h> //For rmdir()
#include <petsc.h>
#include "OutputFileHandler.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestOutputFileHandler : public CxxTest::TestSuite
{
public:

    void TestHandler() throw(Exception)
    {
        // Make a handler that points straight to the CHASTE_TEST_OUTPUT directory.
        OutputFileHandler handler("");
        TS_ASSERT(handler.GetOutputDirectoryFullPath().length() > 0);
        TS_ASSERT_EQUALS(handler.GetOutputDirectoryFullPath(),handler.GetChasteTestOutputDirectory());

        // Make a handler that points straight to a sub-directory.
        std::string dir = "testhandler";
        OutputFileHandler handler2(dir);
        std::string full_dir = handler2.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(full_dir.substr(full_dir.length()-dir.length()-1), dir+"/");
        TS_ASSERT_EQUALS(full_dir, handler2.GetOutputDirectoryFullPath());

        // Check that both can create files
        out_stream p_file_stream;
        p_file_stream = handler.OpenOutputFile("test_file", std::ios::out);
        EXPECT0(system, "test -e " + handler.GetOutputDirectoryFullPath() + "test_file");

        p_file_stream = handler.OpenOutputFile("test_file");
        EXPECT0(system, "test -e " + handler.GetOutputDirectoryFullPath() + "test_file");

        p_file_stream = handler2.OpenOutputFile("test_file");
        EXPECT0(system, "test -e " + handler2.GetOutputDirectoryFullPath() + "test_file");

        p_file_stream = handler2.OpenOutputFile("test_",34,".txt");
        EXPECT0(system, "test -e " + handler2.GetOutputDirectoryFullPath() + "test_34.txt");

        // This should try to write files to /, which isn't allowed (we hope!)
        TS_ASSERT_THROWS_CONTAINS(OutputFileHandler handler3("../../../../../../../../../../../../../../../",false),
                "due to it potentially being above, and cleaning, CHASTE_TEST_OUTPUT.");

        // We don't want other people using CHASTE_TEST_OUTPUT whilst we are messing with it!
        PetscTools::Barrier();

        // Test that the Chaste directory actually influences the location of files
        char *chaste_test_output = getenv("CHASTE_TEST_OUTPUT");

        setenv("CHASTE_TEST_OUTPUT", "", 1/*Overwrite*/);
        OutputFileHandler handler4("whatever");
        rmdir("testoutput/whatever");

        setenv("CHASTE_TEST_OUTPUT", "somewhere_without_trailing_forward_slash", 1/*Overwrite*/);

        OutputFileHandler handler5("whatever");

        rmdir("somewhere_without_trailing_forward_slash/whatever");
        rmdir("somewhere_without_trailing_forward_slash");

        // Reset the location of CHASTE_TEST_OUTPUT
        setenv("CHASTE_TEST_OUTPUT", chaste_test_output, 1/*Overwrite*/);

        // We don't want other people using CHASTE_TEST_OUTPUT whilst we are messing with it!
        PetscTools::Barrier();
    }

    ///todo #1126 test
//    void TestWeCanOnlyDeleteFoldersWeHaveMadeOurselves() throw(Exception)
//    {
//        OutputFileHandler handler("");
//        std::string test_output_path = handler.GetOutputDirectoryFullPath();
//        std::string command;
//
//        if (PetscTools::AmMaster())
//        {
//            command = "mkdir -p " + test_output_path + "cannot_delete_me";
//            system(command.c_str());
//        }
//        // Wait until directory has been created
//        PetscTools::Barrier();
//
//        command = "test -d " + test_output_path + "cannot_delete_me";
//        // Check this folder has been created...
//        TS_ASSERT_EQUALS(system(command.c_str()), 0);
//
//        // Try to use it as an output folder
//        TS_ASSERT_THROWS_CONTAINS(OutputFileHandler bad_handler("cannot_delete_me"),
//                                  "Cannot clean directory");
//
//    }


};

#endif /*TESTOUTPUTFILEHANDLER_HPP_*/
