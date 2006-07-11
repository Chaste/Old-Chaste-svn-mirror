#ifndef TESTINPUTFILEREADER_HPP_
#define TESTINPUTFILEREADER_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>
#include <sstream>
#include "InputFileReader.hpp"

class TestInputFileReader : public CxxTest::TestSuite 
{
public :
    void testInputFileReader(void)
    {
        InputFileReader reader("io/test/data/inputfile.txt");
        
        bool found = false;
        std::string st;
        int n;
        double x,y;
        
        TS_ASSERT_THROWS_ANYTHING(reader.ReadString("", found));

        st = reader.ReadString("MyString", found);
        TS_ASSERT(found==true);
        TS_ASSERT_EQUALS(st,"readmeinstead");

        // won't be read as not in the file - last param says don't quit if not found
        st = reader.ReadString("ShouldntFindThis", found, false);
        TS_ASSERT(found==false);
        TS_ASSERT_EQUALS(st,"");

        // won't be read as not in the commented out
        st = reader.ReadString("AnotherString", found, false);
        TS_ASSERT(found==false);
        TS_ASSERT_EQUALS(st,"");

        // throws because, although 'MyNewString :' is in the file, it is not at
        // the beginning of a line
        TS_ASSERT_THROWS_ANYTHING( reader.ReadString("MyNewString", found) );
        // throws because the NoColon string in the file is not followed by a colon
        TS_ASSERT_THROWS_ANYTHING( reader.ReadString("NoColon", found) );

        n = reader.ReadInt("My Integer", found);
        TS_ASSERT(found==true);
        TS_ASSERT_EQUALS(n,88);

        x = reader.ReadDouble("a double", found);
        TS_ASSERT(found==true);
        TS_ASSERT_DELTA(x,343.435,1e-12);
        
        y = reader.ReadDouble("another double", found);
        TS_ASSERT(found==true);
        TS_ASSERT_DELTA(y,2.4e-4,1e-12);


        // test reading vectors of strings
        std::vector<std::string> strings = reader.ReadVector<std::string>("a vector of strings",2,found);        
        TS_ASSERT(found==true);
        TS_ASSERT_EQUALS(strings.size(),2)
        TS_ASSERT_EQUALS(strings[0],"string1");
        TS_ASSERT_EQUALS(strings[1],"string2");

        // test reading vectors of ints
        std::vector<int> integers = reader.ReadVector<int>("a vector of ints",5,found);        
        TS_ASSERT(found==true);
        TS_ASSERT_EQUALS(integers.size(),5)
        TS_ASSERT_EQUALS(integers[0],3);
        TS_ASSERT_EQUALS(integers[1],4);
        TS_ASSERT_EQUALS(integers[2],5);
        TS_ASSERT_EQUALS(integers[3],2);
        TS_ASSERT_EQUALS(integers[4],77);

        // test reading vectors of doubles
        std::vector<double> doubles = reader.ReadVector<double>("a vector of doubles",4,found);        
        TS_ASSERT(found==true);
        TS_ASSERT_EQUALS(doubles.size(),4)
        TS_ASSERT_DELTA(doubles[0],1,   1e-10);
        TS_ASSERT_DELTA(doubles[1],2.3, 1e-10);
        TS_ASSERT_DELTA(doubles[2],4.2, 1e-10);
        TS_ASSERT_DELTA(doubles[3],1e-3,1e-10);


        // these fields are included more than once in the file, reader should find them twice and
        // throw exceptions
        TS_ASSERT_THROWS_ANYTHING( reader.ReadString("repeated_field_name",found, false) );
        TS_ASSERT_THROWS_ANYTHING( reader.ReadString("another_repeated_field_name",found, false) );


        // ask for too much data, should get exceptions
        TS_ASSERT_THROWS_ANYTHING( reader.ReadVector<std::string>("a vector of strings",2+1,found) );
        TS_ASSERT_THROWS_ANYTHING( reader.ReadVector<int>("a vector of ints",5+1,found) );        
        TS_ASSERT_THROWS_ANYTHING( reader.ReadVector<double>("a vector of doubles",4+1,found) );        

        
        // check nothing is thrown if extra argument saying 'don't quit if not
        // found' is given
        TS_ASSERT_THROWS_NOTHING( reader.ReadString("not there",found,false));
        TS_ASSERT_THROWS_NOTHING( reader.ReadInt   ("not there",found,false));
        TS_ASSERT_THROWS_NOTHING( reader.ReadDouble("not there",found,false));

        TS_ASSERT_THROWS_NOTHING( reader.ReadVector<std::string>("not there",1,found,false));
        TS_ASSERT_THROWS_NOTHING( reader.ReadVector<int>   ("not there",1,found,false));
        TS_ASSERT_THROWS_NOTHING( reader.ReadVector<double>("not there",1,found,false));

        // check nothing is thrown if extra argument saying 'don't quit if not
        // found' is not given
        TS_ASSERT_THROWS_ANYTHING( reader.ReadString("not there",found));
        TS_ASSERT_THROWS_ANYTHING( reader.ReadInt   ("not there",found));
        TS_ASSERT_THROWS_ANYTHING( reader.ReadDouble("not there",found));

        TS_ASSERT_THROWS_ANYTHING( reader.ReadVector<std::string>("not there",1,found));
        TS_ASSERT_THROWS_ANYTHING( reader.ReadVector<int>   ("not there",1,found));
        TS_ASSERT_THROWS_ANYTHING( reader.ReadVector<double>("not there",1,found));

    }
};


#endif /*TESTINPUTFILEREADER_HPP_*/
