#ifndef TESTLOGFILE_HPP_
#define TESTLOGFILE_HPP_

#include <cxxtest/TestSuite.h>
#include "LogFile.hpp"

class TestLogFile : public CxxTest::TestSuite
{
public:
    void TestLogFileCreate()
    {
        LogFile* p_log_file = LogFile::Instance();
    
        // no file set yet
        TS_ASSERT_EQUALS(p_log_file->IsFileSet(), false);
        
        // set the file
        p_log_file->Set(1, "TestLogFile");
        TS_ASSERT_EQUALS(p_log_file->IsFileSet(), true);
        
        // check a new instance works correctly
        LogFile* p_same_log = LogFile::Instance();
        TS_ASSERT_EQUALS(p_same_log->IsFileSet(), true);
    }

    void TestLogStillExists()
    {
        LogFile* p_log = LogFile::Instance();
        TS_ASSERT_EQUALS(p_log->IsFileSet(), true);
    }
    
    void TestClose()
    {
        LogFile::Close();
    
        // check file not set on a new instance
        LogFile* p_log = LogFile::Instance();
        TS_ASSERT_EQUALS(p_log->IsFileSet(), false);
    }

    void TestWritingToFile1()
    {
        LogFile* p_log_file = LogFile::Instance();
        p_log_file->Set(1, "TestLogFile", "log2.txt");
        
        (*p_log_file) << "Some stuff\n" << "Some more\n";
        (*p_log_file) << "Even more\n";
    }

    void TestWritingToFile2()
    {
        LogFile* p_log_file = LogFile::Instance();

        (*p_log_file) << ".. and another bit\n";
        
        (*LogFile::Instance()) << "..and one final bit\n";
        LogFile::Close();
        
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestLogFile/";
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "log2.txt  global/test/data/good_log2.txt").c_str()), 0);
    }
        
    void TestWritingToNoFile()
    {
        LogFile::Close();

        // test no seg faults etc
        (*LogFile::Instance()) << "this won't be written anywhere, as not log file has been created";
    }

    void TestWritingToNewFiles()
    {
        LogFile* p_log_file = LogFile::Instance();
        p_log_file->Set(1, "TestLogFile");
        (*p_log_file) << "data";

        // open a new file without closing the previous
        p_log_file->Set(1, "TestLogFile","log3.txt");
        (*p_log_file) << "data";

        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestLogFile/";
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "log.txt  global/test/data/good_log.txt").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "log3.txt  global/test/data/good_log.txt").c_str()), 0);
        
        LogFile::Close();
    }
    
    
    void TestUsingMacroAndLevels()
    {
        LogFile* p_log_file = LogFile::Instance();
        
        // bad level
        TS_ASSERT_THROWS_ANYTHING(p_log_file->Set( LogFile::MaxLoggingLevel()+1, "TestLogFile") );
        
        p_log_file->Set(1, "TestLogFile", "log4.txt");
        
        unsigned i=0;
        LOG(1, "Level 1 info, will be written. i = " << i); 
        LOG(2, "Level 2 info, WONT be written. i = " << i); 

        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "TestLogFile/";
        
        // this should fail if optimised.
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "log4.txt  global/test/data/good_log4.txt").c_str()), 0);

        LogFile::Close();
    }
};
#endif /*TESTLOGFILE_HPP_*/
