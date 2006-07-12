#ifndef TESTGROUPOFNUMBERSFILEREADER_HPP_
#define TESTGROUPOFNUMBERSFILEREADER_HPP_

#include <cxxtest/TestSuite.h>
#include "GroupOfNumbersFileReader.hpp"

class TestGroupOfNumbersFileReader : public CxxTest::TestSuite 
{
public :
    void testGroupOfNumbersFileReader(void) throw(Exception)
    {
        GroupOfNumbersFileReader<int> reader("io/test/data/some_integer_data");
        std::vector<int> ints = reader.GetData();
        TS_ASSERT_EQUALS(ints.size(), 8);
        TS_ASSERT_EQUALS(ints[0], 0);
        TS_ASSERT_EQUALS(ints[1], 1);
        TS_ASSERT_EQUALS(ints[2], 2);
        TS_ASSERT_EQUALS(ints[3], 3);
        TS_ASSERT_EQUALS(ints[4], 3);
        TS_ASSERT_EQUALS(ints[5], 4);
        TS_ASSERT_EQUALS(ints[6], 5);
        TS_ASSERT_EQUALS(ints[7], 6);

        GroupOfNumbersFileReader<double> d_reader("io/test/data/some_double_data");
        std::vector<double> doubles = d_reader.GetData();
        TS_ASSERT_EQUALS(doubles.size(), 5);
        TS_ASSERT_EQUALS(doubles[0], 0.00001);
        TS_ASSERT_EQUALS(doubles[1], 1.2);
        TS_ASSERT_EQUALS(doubles[2], 2.3);
        TS_ASSERT_EQUALS(doubles[3], 432);
        TS_ASSERT_EQUALS(doubles[4], 1e3);

        // try to read ints from double data - should have error
        TS_ASSERT_THROWS_ANYTHING( GroupOfNumbersFileReader<int> i_reader("io/test/data/some_double_data") );
    }
};


#endif /*TESTGROUPOFNUMBERSFILEREADER_HPP_*/
