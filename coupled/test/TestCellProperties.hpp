#ifndef _TESTCELLPROPERTIES_HPP_
#define _TESTCELLPROPERTIES_HPP_

#include <cxxtest/TestSuite.h>
//#include <iostream>

#include "OdeSolution.hpp"
#include "CellProperties.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"

class TestPhysiologicalProperties : public CxxTest::TestSuite
{
public:
    
    void TestExceptionalBehaviour(void)
    {
        // Check throws an exception if no data given
        std::vector<double> empty;
        TS_ASSERT_THROWS_ANYTHING(CellProperties cell_props(empty, empty));
    }
    
    void TestPhysiologicalPropertiesForRegularLr91(void)
    {
        /*
         * Set stimulus
         */   
        double magnitude_of_stimulus = -80.0;  
        double duration_of_stimulus  = 0.5 ;  // ms                     
        double frequency = 1.0/1000.0; // 1Hz
        double when = 100.0;                                      
        RegularStimulus stimulus(magnitude_of_stimulus,
                                 duration_of_stimulus,
                                 frequency,
                                 when);

        EulerIvpOdeSolver solver;

        /*
         * Solve 
         */
        double start_time = 0.0;   // ms
        double end_time = 3450.0;  // ms
        double time_step = 0.01;   // ms

        LuoRudyIModel1991OdeSystem lr91_ode_system(&solver, time_step, &stimulus);
        
        OdeSolution solution = lr91_ode_system.Compute(start_time, end_time);
        
        // Display solution
//        for (int i=0; i<=solution.GetNumberOfTimeSteps(); i++)
//        {
//            std::cout << solution.mTime[i] << "\t" << solution.mSolutions[i][4]
//                << std::endl;
//        }
        
        // Now calculate the properties
        std::vector<double> voltage=solution.GetVariableAtIndex(4);
        CellProperties  cell_props(voltage, solution.rGetTimes()); // Use default threshold

//        std::cout << "Max upstroke vel: " << cell_props.GetMaxUpstrokeVelocity() << std::endl;
//        std::cout << "Cycle length: " << cell_props.GetCycleLength() << std::endl;
//        std::cout << "Max potential: " << cell_props.GetMaxPotential() << std::endl;
//        std::cout << "Min potential: " << cell_props.GetMinPotential() << std::endl;
//        std::cout << "AP amplitude: " << cell_props.GetActionPotentialAmplitude() << std::endl;
//        std::cout << "APD20: " << cell_props.GetActionPotentialDuration(20) << std::endl;
//        std::cout << "APD50: " << cell_props.GetActionPotentialDuration(50) << std::endl;
//        std::cout << "APD90: " << cell_props.GetActionPotentialDuration(90) << std::endl;
       
        TS_ASSERT_DELTA(cell_props.GetMaxUpstrokeVelocity(), 418.4795, 0.001);
        TS_ASSERT_DELTA(cell_props.GetCycleLength(), 1000.00, 0.01);
        TS_ASSERT_DELTA(cell_props.GetMaxPotential(), 43.1665, 0.0001);
        TS_ASSERT_DELTA(cell_props.GetMinPotential(), -84.4395, 0.0001);
        TS_ASSERT_DELTA(cell_props.GetActionPotentialAmplitude(), 127.606, 0.001);
        TS_ASSERT_DELTA(cell_props.GetActionPotentialDuration(20), 6.66416, 0.00001);
        TS_ASSERT_DELTA(cell_props.GetActionPotentialDuration(50), 271.184, 0.001);
        TS_ASSERT_DELTA(cell_props.GetActionPotentialDuration(90), 361.544, 0.001); // Should use penultimate AP
        TS_ASSERT_DELTA(cell_props.GetTimeAtMaxUpstrokeVelocity(), 3100.7300, 0.001);
    }
};

#endif //_TESTCELLPROPERTIES_HPP_
