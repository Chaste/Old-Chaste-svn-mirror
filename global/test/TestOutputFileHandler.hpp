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

    void TestHandler()
    {
        OutputFileHandler handler("");
        TS_ASSERT(handler.GetOutputDirectoryFullPath("").length() > 0);

        std::string dir = "testhandler";
        OutputFileHandler handler2(dir);
        std::string full_dir = handler2.GetOutputDirectoryFullPath(dir);
        TS_ASSERT_EQUALS(full_dir.substr(full_dir.length()-dir.length()-1), dir+"/");
        TS_ASSERT_EQUALS(full_dir, handler2.GetOutputDirectoryFullPath());

        // Only the master process should write to disk
        if (handler.IsMaster())
        {
            out_stream p_file_stream;
            TS_ASSERT_THROWS_NOTHING(p_file_stream = handler.OpenOutputFile("test_file",
                                                 std::ios::out));

            TS_ASSERT_THROWS_NOTHING(p_file_stream = handler.OpenOutputFile("test_file"));

            TS_ASSERT_THROWS_NOTHING(p_file_stream = handler2.OpenOutputFile("test_file"));

            TS_ASSERT_THROWS_NOTHING(p_file_stream = handler2.OpenOutputFile("test_",34,".txt"));

            // This should try to write files to /, which isn't allowed (we hope!)
            OutputFileHandler handler3("../../../../../../../../../../../../../../../",
                                       false);

            TS_ASSERT_THROWS_THIS(p_file_stream = handler3.OpenOutputFile("test_file"),
                    "Could not open file test_file in " + handler3.GetOutputDirectoryFullPath());
        }

        // Test that the Chaste directory is meaningful, just for coverage purposes

        char *chaste_test_output = getenv("CHASTE_TEST_OUTPUT");

        setenv("CHASTE_TEST_OUTPUT", "", 1/*Overwrite*/);

        handler.GetOutputDirectoryFullPath("whatever");

        rmdir("testoutput/whatever");

        setenv("CHASTE_TEST_OUTPUT", "somewhere_without_trailing_forward_slash", 1/*Overwrite*/);

        handler.GetOutputDirectoryFullPath("whatever");

        rmdir("somewhere_without_trailing_forward_slash/whatever");
        rmdir("somewhere_without_trailing_forward_slash");

        setenv("CHASTE_TEST_OUTPUT", chaste_test_output, 1/*Overwrite*/);
    }

    void TestIsMaster()
    {
        // get an output file handler
        OutputFileHandler handler("");

        PetscTruth is_there;
        PetscInitialized(&is_there);
        TS_ASSERT(is_there);
        PetscInt my_rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);
        if (my_rank==0)
        {
            TS_ASSERT(handler.IsMaster());
        }
        else
        {
            TS_ASSERT(!handler.IsMaster());
        }
    }

};

#endif /*TESTOUTPUTFILEHANDLER_HPP_*/
