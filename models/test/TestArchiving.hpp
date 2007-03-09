#ifndef TESTARCHIVING_HPP_
#define TESTARCHIVING_HPP_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>
#include <cxxtest/TestSuite.h>

#include "OutputFileHandler.hpp"

#include <boost/serialization/vector.hpp>

// see http://www.boost.org/libs/serialization/doc/index.html


class ClassOfSimpleVariables
{
private:    
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>
        archive & mNumber;
        archive & mString;
        archive & mVectorOfDoubles; // include <boost/serialization/vector.hpp> for this
        archive & mVectorOfBools;
    }
    
    int mNumber;
    std::string mString;
    std::vector<double> mVectorOfDoubles;
    std::vector<bool> mVectorOfBools;

public:

    ClassOfSimpleVariables(int initial, 
                           std::string string, 
                           std::vector<double> doubles, 
                           std::vector<bool> bools)
    {
        mNumber = initial;
        mString = string;
        
        mVectorOfDoubles = doubles;
        mVectorOfBools = bools;
    }
    
    int GetNumber() const
    {
        return mNumber;
    }
    
    std::string GetString()
    {
        return mString;
    }
    
    std::vector<double>& GetVectorOfDoubles()
    {
        return mVectorOfDoubles;
    }

    std::vector<bool>& GetVectorOfBools()
    {
        return mVectorOfBools;
    }
};

class TestArchiving : public CxxTest::TestSuite
{
public:
    void TestArchiveSimpleVars()
    {
        OutputFileHandler handler("archive");
        std::string archive_filename;
        archive_filename = handler.GetTestOutputDirectory() + "simple_vars.arch";
        
        // Create an ouput archive 
        {
            std::ofstream ofs(archive_filename.c_str());       
            boost::archive::text_oarchive output_arch(ofs);
                        
            std::vector<double> doubles(3);
            doubles[0] = 1.1;
            doubles[1] = 1.2;
            doubles[2] = 1.3;
                        
            std::vector<bool> bools(2);
            bools[0] = true;
            bools[1] = true;

            ClassOfSimpleVariables i(42,"hello",doubles,bools);

            // cast to const.
            output_arch << static_cast<const ClassOfSimpleVariables&>(i);
        }
        
        {  
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);       
            boost::archive::text_iarchive input_arch(ifs);

            std::vector<double> bad_doubles(1);
            bad_doubles[0] = 10.3;
                        
            std::vector<bool> bad_bools(1);
            bad_bools[0] = false;
    
            ClassOfSimpleVariables j(0,"bye",bad_doubles,bad_bools);

            // read the archive
            input_arch >> j ;

            // Check that the values         
            TS_ASSERT_EQUALS(j.GetNumber(),42);
            TS_ASSERT_EQUALS(j.GetString(),"hello");
            TS_ASSERT_EQUALS(j.GetVectorOfDoubles().size(),3u);
            TS_ASSERT_EQUALS(j.GetVectorOfBools().size(),2u);

            TS_ASSERT_DELTA(j.GetVectorOfDoubles()[0],1.1,1e-12);
            TS_ASSERT_DELTA(j.GetVectorOfDoubles()[1],1.2,1e-12);
            TS_ASSERT_DELTA(j.GetVectorOfDoubles()[2],1.3,1e-12);

            TS_ASSERT_EQUALS(j.GetVectorOfBools()[0],true);
            TS_ASSERT_EQUALS(j.GetVectorOfBools()[1],true);
        }
    }
};


#endif /*TESTARCHIVING_HPP_*/
