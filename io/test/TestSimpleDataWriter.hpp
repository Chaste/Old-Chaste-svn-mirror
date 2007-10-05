#ifndef TESTSIMPLEDATAWRITER_HPP_
#define TESTSIMPLEDATAWRITER_HPP_

#include <cxxtest/TestSuite.h>
#include "SimpleDataWriter.hpp"

class TestSimpleDataWriter : public CxxTest::TestSuite
{
public:
    void TestExceptions()
    {
        std::vector<std::vector<double> > empty_data;
        
        // no data
        TS_ASSERT_THROWS_ANYTHING(SimpleDataWriter bad_writer("SimpleDataWriter", "bad1", empty_data));
        
        std::vector<double> t;
        std::vector<double> x;
        
        t.push_back(0);
        x.push_back(0.1);
        x.push_back(0.2);

        // t and x have different sizes        
        TS_ASSERT_THROWS_ANYTHING(SimpleDataWriter bad_writer("SimpleDataWriter", "bad2", t,x));
    }


    void TestSimpleDataWriterWithStdVecs()
    {
        std::vector<double> t;
        std::vector<double> x;
        std::vector<double> y;
        
        for(unsigned i=0; i<4; i++)
        {
            t.push_back(i);
            x.push_back(2*i);
            y.push_back(i+10);
        }
                
        SimpleDataWriter writer1("SimpleDataWriter", "std_vecs1.dat", t, x);
        
        std::vector<std::vector<double> > full_data;
        full_data.push_back(t);
        full_data.push_back(x);
        full_data.push_back(y);
        
        SimpleDataWriter writer2("SimpleDataWriter", "std_vecs2.dat", full_data, false);
    
        SimpleDataWriter writer3("SimpleDataWriter", "std_vecs3.dat", t, false);
    
        // do the testing now so that to also check the directory wasn't cleaned
        // in the second write       
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "SimpleDataWriter/";
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "std_vecs1.dat  io/test/data/good_std_vec1.dat").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "std_vecs2.dat  io/test/data/good_std_vec2.dat").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "std_vecs3.dat  io/test/data/good_std_vec3.dat").c_str()), 0);
    }
};
#endif /*TESTSIMPLEDATAWRITER_HPP_*/
