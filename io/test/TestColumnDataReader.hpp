#ifndef _TESTCOLUMNDATAREADER_HPP_
#define _TESTCOLUMNDATAREADER_HPP_
// MyTestSuite.h
#include <cxxtest/TestSuite.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "ColumnDataReader.hpp"
#include "global/src/Exception.hpp"

using namespace std;
          
class TestColumnDataReader : public CxxTest::TestSuite 
{

private: 
	ColumnDataReader *mpTestReader;
	
	bool filesMatch(std::string testfileName, std::string goodfileName)
	{	
		bool matching = true;
			
		ifstream testfile(testfileName.c_str(),ios::in);
        ifstream goodfile(goodfileName.c_str(),ios::in);
        std::string teststring;
        std::string goodstring;
        
        if (!testfile.is_open() || !goodfile.is_open())
        {
        	throw new Exception("Files not present.");
        }
        
        while(getline(testfile, teststring))
        {
              getline(goodfile,goodstring);
              if (teststring != goodstring)
              {
              		matching = false;
              }
        }
        
        if(getline(goodfile,goodstring))
        {
        	matching = false;
        }
        
        testfile.close();
        goodfile.close();
        return matching;
	}
	
public:

    void setUp()
    {
		//TS_TRACE("Beginning test...");
    }
    void tearDown()
    {
		// TS_TRACE("Completed test");
    }
    void testCreateColumnReader(void)
    {
        //create a new csvdata writer
        TS_ASSERT_THROWS_ANYTHING(mpTestReader = new ColumnDataReader("testoutput","testdoesnotexist"));
        TS_ASSERT_THROWS_ANYTHING(mpTestReader = new ColumnDataReader("io/test/data","testbad"));
        TS_ASSERT_THROWS_NOTHING(mpTestReader = new ColumnDataReader("testoutput","testfixed"));
        delete mpTestReader; 
    }
};
#endif //_TESTCOLUMNDATAREADER_HPP_
