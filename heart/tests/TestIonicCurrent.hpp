#ifndef _TESTIONICCURRENT_HPP_
#define _TESTIONICCURRENT_HPP_

#include "IonicCurrent.hpp"

class TestIonicCurrent : public CxxTest::TestSuite
{
	public:
	
	// Test ionic current
	void testIonicCurrent(void)
	{
		IonicCurrent *myIonicCurrent;
		
		myIonicCurrent = new IonicCurrent(100.0);
		
		TS_ASSERT(myIonicCurrent->GetMagnitudeOfCurrent() == 100.0);
		
		myIonicCurrent->SetMagnitudeOfCurrent(50.0);
		
		TS_ASSERT(myIonicCurrent->GetMagnitudeOfCurrent() == 50.0);
	}
};

#endif //_TESTIONICCURRENT_HPP_
