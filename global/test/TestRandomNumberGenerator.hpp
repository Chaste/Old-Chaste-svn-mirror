/*

Copyright (C) University of Oxford, 2008

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


#ifndef TESTRANDOMNUMBERS_HPP_
#define TESTRANDOMNUMBERS_HPP_
#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "OutputFileHandler.hpp"
#include "RandomNumberGenerator.hpp"

class TestRandomNumbers : public CxxTest::TestSuite
{
public:
    double ran1;
    
    void TestRandomNumers()
    {
        srandom(0);
        ran1=(double)random()/RAND_MAX;
        
        RandomNumberGenerator *p_gen=RandomNumberGenerator::Instance();
        
        double ran2=p_gen->ranf();
        TS_ASSERT_DELTA(ran1,ran2,1e-7);
        
        RandomNumberGenerator::Destroy();
    }
    
    void TestNewMethodSeed()
    {
        RandomNumberGenerator *p_gen=RandomNumberGenerator::Instance();
        double ran2=p_gen->ranf();
        TS_ASSERT_DELTA(ran1,ran2,1e-7);
        
        RandomNumberGenerator::Destroy();
        
    }
    
    void TestDifferentRandomSeed()
    {
        srandom(36);
        ran1=(double)random()/RAND_MAX;
        
        RandomNumberGenerator *p_gen=RandomNumberGenerator::Instance();
        p_gen->Reseed(36);
        
        
        double ran2=p_gen->ranf();
        TS_ASSERT_DELTA(ran1,ran2,1e-7);
        
        RandomNumberGenerator::Destroy();
        
    }
    
    void TestArchiveRandomNumberGenerator()
    {
        OutputFileHandler handler("archive",false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "random_number.arch";
        
        std::vector<double> generated_numbers;
        
        // Create and archive random number generator
        {	// Save random number generator
        	RandomNumberGenerator *p_gen=RandomNumberGenerator::Instance();
        	p_gen->Reseed(5);

			std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            
            for(unsigned i=0; i<5; i++)
        	{
        		p_gen->ranf();
        	}
        	
        	//A few extra calls before archiving
        	p_gen->ranf();
        	p_gen->randMod(3);        	
        	p_gen->StandardNormalRandomDeviate();
            p_gen->NormalRandomDeviate(0.5, 0.1);
            
            output_arch << static_cast<const RandomNumberGenerator&>(*p_gen);
            
            // Generator saved here - record the next 10 numbers
            
            for(unsigned i=0; i<10; i++)
        	{
        		double random = p_gen->ranf();
        		generated_numbers.push_back(random);
        	}
            
            RandomNumberGenerator::Destroy();
        }
        
        // Restore
        {
            RandomNumberGenerator *p_gen=RandomNumberGenerator::Instance();
        	p_gen->Reseed(25);	// any old seed.
        	for(unsigned i=0; i<7; i++)	// generate some numbers
        	{
        		p_gen->ranf();
        	}
            
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> *p_gen;
            
			// Random Number generator restored. 
			// check it generates the same numbers as the one we saved.
            
			for(unsigned i=0; i<generated_numbers.size(); i++)
        	{
        		double random = p_gen->ranf();
        		TS_ASSERT_DELTA(random,generated_numbers[i],1e-7);
        	}
			
            RandomNumberGenerator::Destroy();
        }
    }
    
    
};

#endif /*TESTRANDOMNUMBERS_HPP_*/
