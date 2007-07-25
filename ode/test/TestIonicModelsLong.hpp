#ifndef _TESTIONICMODELSLONG_HPP_
#define _TESTIONICMODELSLONG_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <ctime>

#include "RunAndCheckIonicModels.hpp"

#include "InitialStimulus.hpp"
#include "RegularStimulus.hpp"

#include "EulerIvpOdeSolver.hpp"

#include "FoxModel2002Modified.hpp"
#include "BackwardEulerFoxModel2002Modified.hpp"

// Note: RunOdeSolverWithIonicModel(), CheckCellModelResults(), CompareCellModelResults()
// are defined in RunAndCheckIonicModels.hpp

class TestIonicModelsLong : public CxxTest::TestSuite
{
public:
    void TestOdeSolverForFox2002WithRegularStimulus(void) throw (Exception)
    {
        clock_t ck_start, ck_end;
        // Set stimulus
        double magnitude = -80.0;
        double duration  = 1.0  ;  // ms
        double start = 50.0; // ms
        double frequency = 1.0/500; // ms^-1
        RegularStimulus stimulus(magnitude, duration, frequency, start);
        
        double end_time = 1000.0; //One second in milliseconds
        double time_step = 0.002;  //2e-6 seconds in milliseconds
                            // 0.005 leads to NaNs.
        
        EulerIvpOdeSolver solver;
        FoxModel2002Modified fox_ode_system(&solver, time_step, &stimulus);
        
        // Solve and write to file
        ck_start = clock();
        RunOdeSolverWithIonicModel(&fox_ode_system,
                                   end_time,
                                   "FoxRegularStimLong",
                                   500);
        ck_end = clock();
        double forward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
                                   
        CheckCellModelResults("FoxRegularStimLong");
        
        // Solve using Backward Euler
        BackwardEulerFoxModel2002Modified backward_system(time_step*5, &stimulus);
        ck_start = clock();
        RunOdeSolverWithIonicModel(&backward_system,
                                   end_time,
                                   "BackwardFoxRegularStimLong",
                                   100);
        ck_end = clock();
        double backward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        
        CompareCellModelResults("FoxRegularStimLong", "BackwardFoxRegularStimLong", 0.15);
        // Mainly for coverage, and to test consistency of GetIIonic
        TS_ASSERT_DELTA(fox_ode_system.GetIIonic(),
                        backward_system.GetIIonic(),
                        1e-6);
        
        std::cout << "Run times:\n\tForward: " << forward
                  << "\n\tBackward: " << backward
                  << std::endl;
        
    }
};


#endif //_TESTIONICMODELSLONG_HPP_
