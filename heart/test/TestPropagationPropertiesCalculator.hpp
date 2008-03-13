#ifndef _TESTPROPAGATIONPROPERTIESCALCULATOR_HPP_
#define _TESTPROPAGATIONPROPERTIESCALCULATOR_HPP_

#include <cxxtest/TestSuite.h>

#include <iostream>

#include "PropagationPropertiesCalculator.hpp"


class TestPropagationPropertiesCalculator : public CxxTest::TestSuite
{
public:

    void TestConductionVelocity1D(void) throw (Exception)
    {
        Hdf5DataReader simulation_data("heart/test/data/MonoDg01d",
                                         "NewMonodomainLR91_1d", false);
        PropagationPropertiesCalculator ppc(&simulation_data);
        double velocity=ppc.CalculateConductionVelocity(5,15,0.1);
        TS_ASSERT_DELTA(velocity, 0.0499, 0.001);
        
        // Should throw because AP does not propagate far enough in simulation time
        TS_ASSERT_THROWS_ANYTHING(ppc.CalculateConductionVelocity(5,95,0.9));
        
        // Should throw because AP does not propagate far enough in simulation time
        TS_ASSERT_THROWS_ANYTHING(ppc.CalculateConductionVelocity(90,100,0.1));
        
        TS_ASSERT_DELTA(ppc.CalculateMaximumUpstrokeVelocity(1),343.9429,0.001);
        
        TS_ASSERT_DELTA(ppc.CalculatePeakMembranePotential(5),23.4467,0.001);
        
        TS_ASSERT_DELTA(ppc.CalculateActionPotentialDuration(50,5),0,0.001);
    }
};

#endif //_TESTPROPAGATIONPROPERTIESCALCULATOR_HPP_
