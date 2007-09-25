#ifndef TESTWNTGRADIENT_HPP_
#define TESTWNTGRADIENT_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include "OutputFileHandler.hpp"
#include "CancerParameters.hpp"
#include "WntGradient.hpp"
#include "SingletonWntGradient.hpp"
#include "WntGradientTypes.hpp"

class TestWntGradient : public CxxTest::TestSuite
{
public:
    void TestWntGradientSetup() throw(Exception)
    {
        CancerParameters::Instance()->Reset();

        WntGradientType this_type = LINEAR;
        
        TS_ASSERT_THROWS_NOTHING(WntGradient wnt_gradient1);
        TS_ASSERT_THROWS_NOTHING(WntGradient wnt_gradient2(this_type));
    }
    
    void TestNoWntGradient() throw(Exception)
    {
        CancerParameters::Instance()->Reset();

        WntGradientType this_type = NONE;
        WntGradient wnt_gradient3(this_type);
        
        double height = 5;
        double wnt_level = 0.0;
        wnt_level = wnt_gradient3.GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);
    }
    
    void TestLinearWntGradient() throw(Exception)
    {
        CancerParameters::Instance()->Reset();

        WntGradientType this_type = LINEAR;
        CancerParameters *params = CancerParameters::Instance();
        WntGradient wnt_gradient3(this_type);
        
        double height = 100;
        double wnt_level = 0.0;
        wnt_level = wnt_gradient3.GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);
        
        height = -1e-12;    // for cells very close to 0 on negative side.
        wnt_level = wnt_gradient3.GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 1.0, 1e-9);        
        
        height = 21.0;
        wnt_level = wnt_gradient3.GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level, 1.0-height/params->GetCryptLength(), 1e-9);
        
        params->SetCryptLength(10.0);
        wnt_level = wnt_gradient3.GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level , 0.0 , 1e-9);
    }
    
    
    void TestOffsetLinearWntGradient() throw(Exception)
    {
        CancerParameters::Instance()->Reset();

        WntGradientType this_type = OFFSET_LINEAR;
        CancerParameters *params = CancerParameters::Instance();
        WntGradient wnt_gradient3(this_type);
        
        double height = 100;
        double wnt_level = 0.0;
        wnt_level = wnt_gradient3.GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);
        
        height = -1e-12;
        wnt_level = wnt_gradient3.GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 1.0, 1e-9);        
        
        height = 21.0;
        wnt_level = wnt_gradient3.GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level, 0.0 , 1e-9);
        
        params->SetCryptLength(10.0);
        wnt_level = wnt_gradient3.GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level , 0.0 , 1e-9);
        // under a third of the way up the crypt.
        params->SetCryptLength(22.0);
        height = 7.0;
        wnt_level = wnt_gradient3.GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level , 1.0 - height/((1.0/3.0)*params->GetCryptLength()) , 1e-9);
        // more than a third of the way up the crypt.
        height = 10.0;
        wnt_level = wnt_gradient3.GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);
    }
    
    void TestArchiveWntGradient()
    {
        CancerParameters::Instance()->Reset();

        OutputFileHandler handler("archive",false);
        std::string archive_filename;
        archive_filename = handler.GetTestOutputDirectory() + "wnt_grad.arch";
        
        // Create an ouput archive
        {
            WntGradientType this_type = LINEAR;
            
            WntGradient* const p_wnt_gradient = new WntGradient(this_type);
            
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            
            output_arch << static_cast<const CancerParameters&>(*CancerParameters::Instance());
            output_arch << p_wnt_gradient;
            
            CancerParameters *inst1 = CancerParameters::Instance();
            TS_ASSERT_DELTA(inst1->GetSG2MDuration(),10.0,1e-12);
            
            delete p_wnt_gradient;
        }
        
        {
            WntGradient* p_wnt;
            
            CancerParameters *inst1 = CancerParameters::Instance();
            
            inst1->SetSG2MDuration(101.0);
            
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            
            // restore from the archive
            input_arch >> *inst1;
            input_arch >> p_wnt;
            
            CancerParameters *inst2 = CancerParameters::Instance();
            TS_ASSERT_EQUALS(inst1, inst2);
            
            // Check
            TS_ASSERT_DELTA(inst1->GetSG2MDuration(),10.0,1e-12);
            double height = 21.0;
            double wnt_level = p_wnt->GetWntLevel(height);
            
            TS_ASSERT_DELTA(wnt_level, 1.0-height/inst1->GetCryptLength(), 1e-9);
        }
    }
    
    
    void TestSingletonWntGradient()
    {
        CancerParameters::Instance()->Reset();
        CancerParameters *params = CancerParameters::Instance();
        
        SingletonWntGradient* p_wnt_gradient = SingletonWntGradient::Instance();
        p_wnt_gradient->SetType(NONE);
        
        double height = 5;
        double wnt_level = 0.0;
        wnt_level = p_wnt_gradient->GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);
        
        TS_ASSERT_THROWS_ANYTHING(p_wnt_gradient->SetType(NONE));
        
        SingletonWntGradient::Destroy();   
 
        p_wnt_gradient = SingletonWntGradient::Instance();

        p_wnt_gradient->SetType(LINEAR);
        
        height = 100;
        wnt_level = 0.0;
        wnt_level = p_wnt_gradient->GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);
        
        height = -1e-12;    // for cells very close to 0 on negative side.
        wnt_level = p_wnt_gradient->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 1.0, 1e-9);        
        
        height = 21.0;
        wnt_level = p_wnt_gradient->GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level, 1.0-height/params->GetCryptLength(), 1e-9);
        
        params->SetCryptLength(10.0);
        wnt_level = p_wnt_gradient->GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level , 0.0 , 1e-9);
    }
};

#endif /*TESTWNTGRADIENT_HPP_*/
