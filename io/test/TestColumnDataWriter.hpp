// MyTestSuite.h
#include <cxxtest/TestSuite.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "ColumnDataWriter.hpp"
#include "global/src/Exception.hpp"

using namespace std;
          
class TestColumnDataWriter : public CxxTest::TestSuite 
{

private: 
    ColumnDataWriter *mpTestWriter;

public:

    void setUp()
    {
        //TS_TRACE("Beginning test...");
    }
    void tearDown()
    {
       // TS_TRACE("Completed test");
    }
    void testCreateColumnWriter(void)
    {
        //create a new csvdata writer
        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("testoutput","test"));
        //check that the output directory exists
        //use the Boost libraries for this check 
    }
    
    void testDefineUnlimitedDimension( void )
    {
        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("testoutput","test"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineUnlimitedDimension("Time","msecs"));
    }
 
    void testDefineFixedDimension( void )
    {
        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("testoutput","test"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineFixedDimension("Node","dimensionless", 5000));
    }
    
    void testDefineVariable( void )
    {
        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("testoutput","test"));
        int ina_var_id = 0;
        int ik_var_id = 0;

        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineUnlimitedDimension("Time","msecs"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineUnlimitedDimension("Time","msecs"));

        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na","milli amperes"));
        TS_ASSERT_THROWS_NOTHING(ik_var_id = mpTestWriter->DefineVariable("I_K","milli amperes"));
        
        TS_ASSERT_EQUALS(ina_var_id, 0);
        TS_ASSERT_EQUALS(ik_var_id, 1);
 
    }

    void testEndDefineMode( void )
    {
        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("testoutput","testdefine"));
        //ending define mode without having defined a dimension and a variable should raise an exception
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->EndDefineMode());
        
        int ina_var_id = 0;
        int ik_var_id = 0;

        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineUnlimitedDimension("Time","msecs"));
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->EndDefineMode());

        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na","milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ik_var_id = mpTestWriter->DefineVariable("I_K","milliamperes"));

        TS_ASSERT_THROWS_NOTHING(mpTestWriter->EndDefineMode());
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->DefineVariable("I_Ca","milli amperes")); 
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->DefineUnlimitedDimension("Time","msecs"));
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->DefineFixedDimension("Node","dimensionless", 5000));
        //mpTestWriter->Close();
    }

    void testPutVariableInUnlimitedFile( void )
    {
        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("testoutput","testunlimited"));
        int time_var_id = 0;
        int ina_var_id = 0;
        int ik_var_id = 0;
        int ica_var_id = 0;
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineUnlimitedDimension("Time","msecs"));
        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na","milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ik_var_id = mpTestWriter->DefineVariable("I_K","milliamperes"));
        TS_ASSERT_THROWS_NOTHING(time_var_id = mpTestWriter->DefineVariable("Time","msecs"));
        TS_ASSERT_THROWS_NOTHING(ica_var_id = mpTestWriter->DefineVariable("I_Ca","milliamperes"));
        //TS_TRACE("Here");
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->EndDefineMode());
        int i = 12;

        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(time_var_id, 0.1));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ina_var_id, (double) i));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ica_var_id, 33.124));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ik_var_id, 7124.12355553));
        mpTestWriter->AdvanceAlongUnlimitedDimension();
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ica_var_id, 63.124));
        mpTestWriter->Close();
        //TODO: compare output with testcompare.dat
        
    }
    
    
    // Note to Dan: Gary wrote this one, it might not work.
    void testPutVariableInUnlimitedNegativeFile( void )
    {
        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("testoutput","testunlimitednegative"));
        int time_var_id = 0;
        int ina_var_id = 0;
        int ik_var_id = 0;
        int ica_var_id = 0;
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineUnlimitedDimension("Time","msecs"));
        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na","milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ik_var_id = mpTestWriter->DefineVariable("I_K","milliamperes"));
        TS_ASSERT_THROWS_NOTHING(time_var_id = mpTestWriter->DefineVariable("Time","msecs"));
        TS_ASSERT_THROWS_NOTHING(ica_var_id = mpTestWriter->DefineVariable("I_Ca","milliamperes"));
        //TS_TRACE("Here");
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->EndDefineMode());
        int i = 12;

        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(time_var_id, 0.2));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ina_var_id, (double) -i));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ica_var_id, 33.124));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ik_var_id, 7124.12355553));
        mpTestWriter->AdvanceAlongUnlimitedDimension();
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ica_var_id, -33.124));
        mpTestWriter->Close();
        
    }

    void testPutVariableInFixedFile( void )
    {

        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("testoutput","testfixed"));
        int node_var_id = 0;
        int ina_var_id = 0;
        int ik_var_id = 0;
        int ica_var_id = 0;
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineFixedDimension("Node","dimensionless",4));
        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na","milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ik_var_id = mpTestWriter->DefineVariable("I_K","milliamperes"));
        TS_ASSERT_THROWS_NOTHING(node_var_id = mpTestWriter->DefineVariable("Node","dimensionless")) ;
        TS_ASSERT_THROWS_NOTHING(ica_var_id = mpTestWriter->DefineVariable("I_Ca","milliamperes"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->EndDefineMode());
        int i = 12;

        mpTestWriter->PutVariable(ina_var_id, (double) i,0);
        mpTestWriter->PutVariable(ina_var_id, (double) i,1);
        mpTestWriter->PutVariable(ica_var_id, 33.124,3);
        mpTestWriter->PutVariable(ik_var_id, 7124.12355553,3);
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->AdvanceAlongUnlimitedDimension());
        mpTestWriter->PutVariable(ica_var_id, 63.124,2);
        mpTestWriter->PutVariable(node_var_id, 1,0);
        mpTestWriter->PutVariable(node_var_id, 4,3);
        mpTestWriter->Close();
    }
    
    
    void testPutNegativeVariableInFixedFile( void )
    {

        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("testoutput","testfixed_negatives"));
        int node_var_id = 0;
        int ina_var_id = 0;
        int ik_var_id = 0;
        int ica_var_id = 0;
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineFixedDimension("Node","dimensionless",4));
        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na","milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ik_var_id = mpTestWriter->DefineVariable("I_K","milliamperes"));
        TS_ASSERT_THROWS_NOTHING(node_var_id = mpTestWriter->DefineVariable("Node","dimensionless")) ;
        TS_ASSERT_THROWS_NOTHING(ica_var_id = mpTestWriter->DefineVariable("I_Ca","milliamperes"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->EndDefineMode());
        int i = 12;

        mpTestWriter->PutVariable(ina_var_id, (double) i,0);
        mpTestWriter->PutVariable(ina_var_id, (double) -i,1);
        mpTestWriter->PutVariable(ica_var_id, -33.124,3);
        mpTestWriter->PutVariable(ik_var_id, 7124.12355553,3);
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->AdvanceAlongUnlimitedDimension());
        mpTestWriter->PutVariable(ica_var_id, -63.124,2);
        mpTestWriter->PutVariable(node_var_id, 1,0);
        mpTestWriter->PutVariable(node_var_id, -4,3);
        mpTestWriter->Close();
    }
    
    void testPutVariableInFixedandUnlimitedFile( void )
    {

        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("testoutput","testfixedandunlimited"));
        int time_var_id = 0;
        int ina_var_id = 0;
        int ik_var_id = 0;
        int ica_var_id = 0;
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineFixedDimension("Node","dimensionless",4));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineUnlimitedDimension("Time","msecs"));
        TS_ASSERT_THROWS_NOTHING(time_var_id = mpTestWriter->DefineVariable("Time","msecs"));
        TS_ASSERT_THROWS_NOTHING(ina_var_id = mpTestWriter->DefineVariable("I_Na","milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ik_var_id = mpTestWriter->DefineVariable("I_K","milliamperes"));
        TS_ASSERT_THROWS_NOTHING(ica_var_id = mpTestWriter->DefineVariable("I_Ca","milliamperes"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->EndDefineMode());
        int i = 12;

        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(time_var_id, 0.1));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ina_var_id, (double) i,0));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ina_var_id, (double) i,1));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ica_var_id, 33.124,3));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ik_var_id, 7124.12355553,3));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->AdvanceAlongUnlimitedDimension());
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(ica_var_id, 63.124,2));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->PutVariable(time_var_id, 0.2));
        TS_ASSERT_THROWS_ANYTHING(mpTestWriter->PutVariable(time_var_id, 0.2,3));
        mpTestWriter->Close();
    }
    
    void testDefineFilesMatch( void )
    {
        ifstream testfile("testoutput/testdefine.dat",ios::in);
        ifstream goodfile("io/test/data/testdefine_good.dat",ios::in);
        std::string teststring;
        std::string goodstring;
        while(getline(testfile, teststring))
        {
              getline(goodfile,goodstring);
              TS_ASSERT_EQUALS(teststring,goodstring);
        }
        testfile.close();
        goodfile.close();
        
    }
    void testFixedFilesMatch( void )
    {
        ifstream testfile("testoutput/testfixed.dat",ios::in);
        ifstream goodfile("io/test/data/testfixed_good.dat",ios::in);
        std::string teststring;
        std::string goodstring;
        while(getline(goodfile, goodstring))
        {
              getline(testfile,teststring);
              TS_ASSERT_EQUALS(teststring,goodstring);
        }
        testfile.close();
        goodfile.close();
        
    }
    void testFixedAndUnlimitedTimeFilesMatch( void )
    {
        ifstream testfile("testoutput/testfixedandunlimitedTime.dat",ios::in);
        ifstream goodfile("io/test/data/testfixedandunlimitedTime_good.dat",ios::in);
        std::string teststring;
        std::string goodstring;
        while(getline(goodfile, goodstring))
        {
              getline(testfile,teststring);
              TS_ASSERT_EQUALS(teststring,goodstring);
        }
        testfile.close();
        goodfile.close();
        
    }
    void testFixedAndUnlimitedFilesMatch( void )
    {
        ifstream testfile("testoutput/testfixedandunlimited.dat",ios::in);
        ifstream goodfile("io/test/data/testfixedandunlimited_good.dat",ios::in);
        std::string teststring;
        std::string goodstring;
        while(getline(goodfile, goodstring))
        {
              getline(testfile,teststring);
              TS_ASSERT_EQUALS(teststring,goodstring);
        }
        testfile.close();
        goodfile.close();
        
    }
    void testUnlimitedFilesMatch( void )
    {
        ifstream testfile("testoutput/testunlimited.dat",ios::in);
        ifstream goodfile("io/test/data/testunlimited_good.dat",ios::in);
        std::string teststring;
        std::string goodstring;
        while(getline(goodfile, goodstring))
        {
              getline(testfile,teststring);
              TS_ASSERT_EQUALS(teststring,goodstring);
        }
        testfile.close();
        goodfile.close();
        
    }
    
    //Gary wrote these two
    void testUnlimitedNegativeFilesMatch( void )
    {
        ifstream testfile("testoutput/testunlimitednegative.dat",ios::in);
        ifstream goodfile("io/test/data/testunlimitednegative_good.dat",ios::in);
        std::string teststring;
        std::string goodstring;
        while(getline(goodfile, goodstring))
        {
              getline(testfile,teststring);
              TS_ASSERT_EQUALS(teststring,goodstring);
        }
        testfile.close();
        goodfile.close();
	}
	
	void testFixedNegativeFilesMatch( void )
    {
        ifstream testfile("testoutput/testfixed_negatives.dat",ios::in);
        ifstream goodfile("io/test/testfixed_negatives_good.dat",ios::in);
        std::string teststring;
        std::string goodstring;
        while(getline(goodfile, goodstring))
        {
              getline(testfile,teststring);
              TS_ASSERT_EQUALS(teststring,goodstring);
        }
        testfile.close();
        goodfile.close();
	}
    
};
