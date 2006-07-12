#ifndef TESTOUTPUTFILEHANDLER_HPP_
#define TESTOUTPUTFILEHANDLER_HPP_

#include <cxxtest/TestSuite.h>
#include <string>
#include <fstream>
#include "OutputFileHandler.hpp"

class TestOutputFileHandler : public CxxTest::TestSuite 
{
public:

	void TestHandler(void)
	{
		OutputFileHandler handler("");
		TS_ASSERT(handler.GetTestOutputDirectory("").length() > 0);
		
		std::string dir = "testhandler";
		OutputFileHandler handler2(dir);
		std::string full_dir = handler2.GetTestOutputDirectory(dir);
		TS_ASSERT_EQUALS(full_dir.substr(full_dir.length()-dir.length()-1), dir+"/");
        TS_ASSERT_EQUALS(full_dir, handler2.GetTestOutputDirectory());
		
		out_stream p_file_stream;
        TS_ASSERT_THROWS_NOTHING(p_file_stream = handler.OpenOutputFile("test_file",
                                                                        std::ios::out));
        
        TS_ASSERT_THROWS_NOTHING(p_file_stream = handler.OpenOutputFile("test_file"));
        
        TS_ASSERT_THROWS_NOTHING(p_file_stream = handler2.OpenOutputFile("test_file"));
        
        // This should try to write files to /, which isn't allowed (we hope!)
        OutputFileHandler handler3("../../../../../../../../../../../../../../../",
                                   false);
        TS_ASSERT_THROWS_ANYTHING(p_file_stream = handler3.OpenOutputFile("test_file"));
    }
};

#endif /*TESTOUTPUTFILEHANDLER_HPP_*/
