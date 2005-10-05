#ifndef _TESTPHYSIOLOGICALPROPERTIES_HPP_
#define _TESTPHYSIOLOGICALPROPERTIES_HPP_

#include <cxxtest/TestSuite.h>

#include <iostream>

#include "OdeSolution.hpp"
#include "PhysiologicalProperties.hpp"

#include "RegularStimulus.hpp"

#include "EulerIvpOdeSolver.hpp"

#include "AbstractOdeSystem.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"

class TestPhysiologicalProperties : public CxxTest::TestSuite
{
public:
    
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

        LuoRudyIModel1991OdeSystem lr91_ode_system(&stimulus);

        EulerIvpOdeSolver solver;
        
        /*
         * Solve 
         */
        double start_time = 0.0;   // ms
        double end_time = 3500.0;   // ms
        double time_step = 0.01;   // ms
                
        OdeSolution solution = solver.Solve(&lr91_ode_system,
                                            start_time,
                                            end_time,
                                            time_step,
                                            lr91_ode_system.mInitialConditions);
        
        // Display solution
//        for (int i=0; i<=solution.GetNumberOfTimeSteps(); i++)
//        {
//            std::cout << solution.mTime[i] << "\t" << solution.mSolutions[i][4]
//                << std::endl;
//        }
        
        // Now calculate the properties
        std::vector<double> voltage=solution.GetVariableAtIndex(4);
        PhysiologicalProperties  pp(voltage, solution.mTime); // Use default threshold

//        std::cout << "Max upstroke vel: " << pp.GetMaxUpstrokeVelocity() << std::endl;
//        std::cout << "Cycle length: " << pp.GetCycleLength() << std::endl;
//        std::cout << "Max potential: " << pp.GetMaxPotential() << std::endl;
//        std::cout << "Min potential: " << pp.GetMinPotential() << std::endl;
//        std::cout << "AP amplitude: " << pp.GetActionPotentialAmplitude() << std::endl;
//        std::cout << "APD20: " << pp.GetActionPotentialDuration(20) << std::endl;
//        std::cout << "APD50: " << pp.GetActionPotentialDuration(50) << std::endl;
//        std::cout << "APD90: " << pp.GetActionPotentialDuration(90) << std::endl;
       
        TS_ASSERT_DELTA(pp.GetMaxUpstrokeVelocity(), 418.834, 0.001);
        TS_ASSERT_DELTA(pp.GetCycleLength(), 999.99, 0.01);
        TS_ASSERT_DELTA(pp.GetMaxPotential(), 43.1772, 0.0001);
        TS_ASSERT_DELTA(pp.GetMinPotential(), -84.4392, 0.0001);
        TS_ASSERT_DELTA(pp.GetActionPotentialAmplitude(), 127.616, 0.001);
        TS_ASSERT_DELTA(pp.GetActionPotentialDuration(20), 6.65597, 0.00001);
        TS_ASSERT_DELTA(pp.GetActionPotentialDuration(50), 271.156, 0.001);
        TS_ASSERT_DELTA(pp.GetActionPotentialDuration(90), 361.546, 0.001);
        TS_ASSERT_DELTA(pp.GetTimeAtMaxUpstrokeVelocity(), 3100.7400, 0.001);
    }
};

#endif //_TESTPHYSIOLOGICALPROPERTIES_HPP_
