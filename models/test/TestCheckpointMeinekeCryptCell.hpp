#ifndef TESTCHECKPOINTMEINEKECRYPTCELL_HPP_
#define TESTCHECKPOINTMEINEKECRYPTCELL_HPP_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>
#include <cxxtest/TestSuite.h>

#include "OutputFileHandler.hpp"

class Number
{
private:    
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & number;
    }

public:
    int number;
    Number(int initial)
    {
        number = initial;
    }
};

class TestCheckpointMeinekeCryptCell : public CxxTest::TestSuite
{
public:
    void TestArchiveInt()
    {

        // Create an ouput archive
        
        OutputFileHandler handler("archive");
        std::string archive_filename;
        archive_filename = handler.GetTestOutputDirectory() + "int.arch";
        
        {
        std::ofstream ofs(archive_filename.c_str());       
        boost::archive::text_oarchive output_arch(ofs);
        
        // write an integer of value 42 to the archive
        const Number i(42);
        output_arch << i;
        }
        
        {  
        // Create an input archive
        std::ifstream ifs(archive_filename.c_str(), std::ios::binary);       
        boost::archive::text_iarchive input_arch(ifs);
        
        // read the archive
        Number j(0);
        input_arch >> j ;
        // Check that the value of 42 was read from the archive        
        TS_ASSERT_EQUALS(j.number,42);
        }
    }
};


#endif /*TESTCHECKPOINTMEINEKECRYPTCELL_HPP_*/
