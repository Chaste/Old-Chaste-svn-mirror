#ifndef TESTOUTPUTFILEHANDLER_HPP_
#define TESTOUTPUTFILEHANDLER_HPP_

#include <cxxtest/TestSuite.h>
#include <string>
#include <fstream>
#include <unistd.h> //For rmdir()
#include <petsc.h>
#include "OutputFileHandler.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestOutputFileHandler : public CxxTest::TestSuite
{
public:

    void TestHandler(void)
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
            TS_ASSERT_THROWS_ANYTHING(p_file_stream = handler3.OpenOutputFile("test_file"));
        }
        
        // Test that the Chaste directory is meaningful, just for coverage purposes
        
        char *chaste_test_output = getenv("CHASTE_TEST_OUTPUT");
        
        setenv("CHASTE_TEST_OUTPUT", "", 1/*Overwrite*/);
        
        handler.GetOutputDirectoryFullPath("whatever");
        
        rmdir("whatever");
        
        setenv("CHASTE_TEST_OUTPUT", "somewhere_without_trailing_forward_slash", 1/*Overwrite*/);
        
        handler.GetOutputDirectoryFullPath("whatever");
        
        rmdir("somewhere_without_trailing_forward_slash/whatever");
        rmdir("somewhere_without_trailing_forward_slash");
        
        setenv("CHASTE_TEST_OUTPUT", chaste_test_output, 1/*Overwrite*/);
    }
    
    void TestIsMaster(void)
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
