#ifndef TESTWNTGRADIENT_HPP_
#define TESTWNTGRADIENT_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include "OutputFileHandler.hpp"
#include "CancerParameters.hpp"
#include "WntGradient.hpp"
#include "WntGradientTypes.hpp"

class TestWntGradient : public CxxTest::TestSuite
{
public:
    void TestWntGradientSetup() throw(Exception)
    {
    	WntGradientType this_type = LINEAR;
        
        TS_ASSERT_THROWS_NOTHING(WntGradient wnt_gradient1);
        
        TS_ASSERT_THROWS_NOTHING(WntGradient wnt_gradient2(this_type));
    }
    
    void TestNoWntGradient() throw(Exception)
    {    
    	WntGradientType this_type = NONE;
        //CancerParameters *params = CancerParameters::Instance();
		WntGradient wnt_gradient3(this_type);

		double height = 5;
		double wnt_level = 0.0;
		wnt_level = wnt_gradient3.GetWntLevel(height);

		TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);
    }
    
    void TestLinearWntGradient() throw(Exception)
    {    
    	WntGradientType this_type = LINEAR;
        CancerParameters *params = CancerParameters::Instance();
		WntGradient wnt_gradient3(this_type);

		double height = 100;
		double wnt_level = 0.0;
		wnt_level = wnt_gradient3.GetWntLevel(height);

		TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);

		height = 21.0;
		wnt_level = wnt_gradient3.GetWntLevel(height);
		
		TS_ASSERT_DELTA(wnt_level, 1.0-height/params->GetCryptLength(), 1e-9);
        
        params->SetCryptLength(10.0);
        wnt_level = wnt_gradient3.GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level , 0.0 , 1e-9);
        
    }
    
    
    void TestOffsetLinearWntGradient() throw(Exception)
    {    
    	WntGradientType this_type = OFFSET_LINEAR;
        CancerParameters *params = CancerParameters::Instance();
		WntGradient wnt_gradient3(this_type);

		double height = 100;
		double wnt_level = 0.0;
		wnt_level = wnt_gradient3.GetWntLevel(height);

		TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);

		height = 21.0;
		wnt_level = wnt_gradient3.GetWntLevel(height);
		
		TS_ASSERT_DELTA(wnt_level, 0.0 , 1e-9);
        
        params->SetCryptLength(10.0);
        wnt_level = wnt_gradient3.GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level , 0.0 , 1e-9);
        
        params->SetCryptLength(22.0);
        height = 10.0;
        wnt_level = wnt_gradient3.GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level , 1.0 - 1.5*height/params->GetCryptLength() , 1e-9);
        
    }
    
    void TestArchiveWntGradient()
    {
        OutputFileHandler handler("archive");
        std::string archive_filename;
        archive_filename = handler.GetTestOutputDirectory() + "wnt_grad.arch";
        
        // Create an ouput archive 
        {
            WntGradientType this_type = LINEAR;
            
            WntGradient wnt_gradient(this_type);
    
            std::ofstream ofs(archive_filename.c_str());       
            boost::archive::text_oarchive output_arch(ofs);
            
            output_arch << static_cast<const CancerParameters&>(*CancerParameters::Instance());
            output_arch << static_cast<const WntGradient&>(wnt_gradient);
            
            CancerParameters *inst1 = CancerParameters::Instance();
            TS_ASSERT_DELTA(inst1->GetSG2MDuration(),10.0,1e-12);
        }
        
        {  
            WntGradientType this_type = NONE;
            WntGradient wnt_gradient(this_type);
            
            CancerParameters *inst1 = CancerParameters::Instance();
            
            inst1->SetSG2MDuration(101.0);
            
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);       
            boost::archive::text_iarchive input_arch(ifs);
            
            // restore from the archive
            input_arch >> *inst1;
            input_arch >> wnt_gradient;
            
            CancerParameters *inst2 = CancerParameters::Instance();
            TS_ASSERT_EQUALS(inst1, inst2);
            
            // Check 
            TS_ASSERT_DELTA(inst1->GetSG2MDuration(),10.0,1e-12);
            double height = 21.0;
            double wnt_level = wnt_gradient.GetWntLevel(height);
        
            TS_ASSERT_DELTA(wnt_level, 1.0-height/inst1->GetCryptLength(), 1e-9);
        }
    }
};

#endif /*TESTWNTGRADIENT_HPP_*/
