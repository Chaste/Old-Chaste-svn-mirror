#ifndef TESTFABERRUDYAGAINSTLIMPET_HPP_
#define TESTFABERRUDYAGAINSTLIMPET_HPP_

#include <cxxtest/TestSuite.h>
#include "RegularStimulus.hpp"
#include "FaberRudy2000Version3Optimised.hpp"
#include "EulerIvpOdeSolver.hpp"

class TestFaberRudyAgainstLimpet : public CxxTest::TestSuite
{
public:
    void TestFaberRudyOptimised(void) throw(Exception)
    {

        
        RegularStimulus stimulus(-50.0,     // magnitude uA/cm^2
                                 1.0,      // duration ms
                                 1000.0,   // period ms
                                 0.0);     // start time ms
        
        EulerIvpOdeSolver solver;
        double time_step = 0.001;                         
        
        FaberRudy2000Version3Optimised cell_model(&solver, time_step, &stimulus);
        
        
        //OdeSolution solution =
        cell_model.ComputeExceptVoltage(0.0,     // start time
                                         100.0); // endTime
                                                  
//        solution.WriteToFile("TestFRvLimpet",
//                             "FaberRudy",
//                             &cell_model,
//                             "ms",
//                             1000u, // steps per row
//                             false);
    }
};

#endif /*TESTFABERRUDYAGAINSTLIMPET_HPP_*/
