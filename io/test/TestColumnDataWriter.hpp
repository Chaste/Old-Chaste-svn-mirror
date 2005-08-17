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
    void testCreateColumnWriter(void)
    {
        //create a new csvdata writer
        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("testoutput","test"));
        //check that the output directory exists
        //use the Boost libraries for this check
        
        delete mpTestWriter; 
    }
    
    void testDefineUnlimitedDimension( void )
    {
        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("testoutput","test"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineUnlimitedDimension("Time","msecs"));
        
        delete mpTestWriter;
    }
 
    void testDefineFixedDimension( void )
    {
        TS_ASSERT_THROWS_NOTHING(mpTestWriter = new ColumnDataWriter("testoutput","test"));
        TS_ASSERT_THROWS_NOTHING(mpTestWriter->DefineFixedDimension("Node","dimensionless", 5000));
        
        delete mpTestWriter;
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

		delete mpTestWriter;
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
        
        delete mpTestWriter;

		TS_ASSERT(filesMatch("testoutput/testdefine.dat", 
		                     "io/test/data/testdefine_good.dat"));
		                     
		TS_ASSERT(filesMatch("testoutput/testdefine.info", 
		                     "io/test/data/testdefine_good.info"));
		                     
		TS_ASSERT(!filesMatch("testoutput/testdefine.info", 
		                      "io/test/data/testdefine_bad.info"));
        
        
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
        
		delete mpTestWriter;

		TS_ASSERT(filesMatch("testoutput/testunlimited.dat", 
		                     "io/test/data/testunlimited_good.dat"));
		                     
    }
    
    
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
    
    	delete mpTestWriter; 

		TS_ASSERT(filesMatch("testoutput/testunlimitednegative.dat", 
		                     "io/test/data/testunlimitednegative_good.dat"));
		                      
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

		delete mpTestWriter;

		TS_ASSERT(filesMatch("testoutput/testfixed.dat", 
		                     "io/test/data/testfixed_good.dat"));
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
    
    	delete mpTestWriter;
    	
		TS_ASSERT(filesMatch("testoutput/testfixed_negatives.dat", 
		                     "io/test/data/testfixed_negatives_good.dat"));

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

		delete mpTestWriter;
		
		TS_ASSERT(filesMatch("testoutput/testfixedandunlimitedTime.dat", 
		                     "io/test/data/testfixedandunlimitedTime_good.dat"));

		TS_ASSERT(filesMatch("testoutput/testfixedandunlimited_1.dat", 
		                     "io/test/data/testfixedandunlimited_1_good.dat"));
    }
    
};
