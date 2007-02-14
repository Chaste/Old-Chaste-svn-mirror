#ifndef TESTWNTGRADIENT_HPP_
#define TESTWNTGRADIENT_HPP_

#include <cxxtest/TestSuite.h>

#include "CancerParameters.hpp"
#include "WntGradient.hpp"
#include "WntGradientTypes.hpp"

class TestWntGradient : public CxxTest::TestSuite
{
public:
    void TestWntGradientSetup() throw(Exception)
    {
    	WntGradientType this_type = LINEAR;
        
        TS_ASSERT_THROWS_ANYTHING(WntGradient wnt_gradient1);
        
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
};

#endif /*TESTWNTGRADIENT_HPP_*/
