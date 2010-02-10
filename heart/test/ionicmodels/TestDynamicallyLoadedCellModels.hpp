/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef TESTDYNAMICALLYLOADEDCELLMODELS_HPP_
#define TESTDYNAMICALLYLOADEDCELLMODELS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/shared_ptr.hpp>

#include "RunAndCheckIonicModels.hpp"
//#include "LuoRudyIModel1991OdeSystem.hpp"
#include "SimpleStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "DynamicCellModelLoader.hpp"
#include "ChasteBuildRoot.hpp"

class TestDynamicallyLoadedCellModels : public CxxTest::TestSuite
{
public:
    /**
     * This is based on TestOdeSolverForLR91WithDelayedSimpleStimulus from
     * TestIonicModels.hpp.
     */
    void TestDynamicallyLoadedLr91(void) throw(Exception)
    {
        // Set stimulus
        double magnitude = -25.5;
        double duration  = 2.0 ;  // ms
        double when = 50.0; // ms

        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude, duration, when));
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        double end_time = 1000.0; //One second in milliseconds

        // Load the cell model dynamically
        std::string model_name = "libDynamicallyLoadableLr91.so";
        DynamicCellModelLoader loader(ChasteComponentBuildDir("heart") + "dynamic/" + model_name);
        AbstractCardiacCell* p_cell = loader.CreateCell(p_solver, p_stimulus);

        // Simple sanity check
        TS_ASSERT_EQUALS(p_cell->GetVoltageIndex(), 4u);

        // Solve and write to file
        clock_t ck_start = clock();
        RunOdeSolverWithIonicModel(p_cell,
                                   end_time,
                                   "DynamicallyLoadableLr91");
        clock_t ck_end = clock();
        double forward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        std::cout << "\n\tSolve time: " << forward << std::endl;

        // Compare with 'normal' LR91 model results
        CheckCellModelResults("DynamicallyLoadableLr91", "Lr91DelayedStim");

        // Test GetIIonic against hardcoded result from TestIonicModels.hpp
        RunOdeSolverWithIonicModel(p_cell,
                                   60.0,
                                   "DynamicallyLoadableLr91GetIIonic");
        TS_ASSERT_DELTA(p_cell->GetIIonic(), 1.9411, 1e-3);

        // Need to delete cell model
        delete p_cell;
    }
    
    void TestExceptions() throw(Exception)
    {
        // Try loading a .so that doesn't exist
        std::string file_name = "non-existent-file-we-hope";
        TS_ASSERT_THROWS_CONTAINS(DynamicCellModelLoader loader(file_name),
                                  "Unable to load .so file '" + file_name + "':");

        // Try loading a .so that doesn't define a cell model
        file_name = "libNotACellModel.so";
        TS_ASSERT_THROWS_CONTAINS(DynamicCellModelLoader loader(ChasteComponentBuildDir("heart") + "dynamic/" + file_name),
                                  "Failed to load cell creation function from .so file");
    }
};

#endif /* TESTDYNAMICALLYLOADEDCELLMODELS_HPP_ */
