#ifndef TESTCANCERPARAMETERS_HPP_
#define TESTCANCERPARAMETERS_HPP_

#include <cxxtest/TestSuite.h>

#include "CancerParameters.hpp"

class TestCancerParameters : public CxxTest::TestSuite
{
public:
    void TestGettersAndSetters()
    {
        CancerParameters *inst1 = CancerParameters::Instance();
        
        inst1->SetSG2MDuration(11.0);
        inst1->SetStemCellCycleTime(-25.0);
        inst1->SetTransitCellCycleTime(-35.0);
        inst1->SetMaxTransitGenerations(666u);
        inst1->SetCryptLength(-1.0);
        inst1->SetMeinekeLambda(-2.0);
        
        CancerParameters *inst2 = CancerParameters::Instance();
        
        TS_ASSERT_DELTA(inst2->GetSG2MDuration(), 11.0 , 1e-12);
        TS_ASSERT_DELTA(inst2->GetStemCellCycleTime(), -25.0, 1e-12);
        TS_ASSERT_DELTA(inst2->GetTransitCellCycleTime(), -35.0, 1e-12);
        TS_ASSERT_EQUALS(inst2->GetMaxTransitGenerations(), 666u);
        TS_ASSERT_DELTA(inst2->GetCryptLength(), -1.0, 1e-12);
        TS_ASSERT_DELTA(inst2->GetMeinekeLambda(), -2.0, 1e-12);
    }
};

#endif /*TESTCANCERPARAMETERS_HPP_*/
