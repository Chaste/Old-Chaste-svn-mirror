#ifndef _TESTPROPAGATIONPROPERTIESCALCULATOR_HPP_
#define _TESTPROPAGATIONPROPERTIESCALCULATOR_HPP_

// Element.hpp includes the Boost ublas objects - these need to
// be included early...  We think.  We're not that sure.
#include "Element.hpp"

#include <cxxtest/TestSuite.h>

#include <iostream>

#include "PropagationPropertiesCalculator.hpp"
#include "ColumnDataReader.hpp"

class TestPropagationPropertiesCalculator : public CxxTest::TestSuite
{
public:

    void TestConductionVelocity1D(void)
    {
        ColumnDataReader simulation_data("coupled/test/data/MonoDg01d",
                                         "NewMonodomainLR91_1d", false);
        PropagationPropertiesCalculator ppc(&simulation_data);
        double velocity=ppc.CalculateConductionVelocity(5,15,0.1);
        TS_ASSERT_DELTA(velocity, 0.0961538, 0.001);
        
        // Should throw because AP does not propagate far enough in simulation time
        TS_ASSERT_THROWS_ANYTHING(ppc.CalculateConductionVelocity(5,95,0.9));
        
        // Should throw because AP does not propagate far enough in simulation time
        TS_ASSERT_THROWS_ANYTHING(ppc.CalculateConductionVelocity(90,100,0.1));
    }
};

#endif //_TESTPROPAGATIONPROPERTIESCALCULATOR_HPP_
