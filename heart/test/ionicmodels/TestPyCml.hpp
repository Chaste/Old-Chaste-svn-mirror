/*

Copyright (C) University of Oxford, 2005-2009

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


#ifndef _TESTPYCML_HPP_
#define _TESTPYCML_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>
#include <string>
#include <ctime>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "RunAndCheckIonicModels.hpp"

#include "Exception.hpp"
#include "SimpleStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "ArchiveLocationInfo.hpp"

#include "luo_rudy_1991.hpp"
#include "luo_rudy_1991Opt.hpp"


class TestPyCml : public CxxTest::TestSuite
{
public:
     /**
      * This test is designed to quickly check that PyCml-generated code matches the Chaste interfaces,
      * and gives expected results.
      * 
      * \todo #1030 run PyCml automatically, rather than having to generate the .hpp files by hand.
      */
     void TestPyCmlCodeGeneration()
     {
        clock_t ck_start, ck_end;

        //
        // Set up cells
        //
        double magnitude = -25.5;
        double duration  = 2.0 ;  // ms
        double when = 50.0; // ms
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude, duration, when));
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        double end_time = 1000.0; //One second in milliseconds

        // Normal model
        Cellluo_rudy_1991FromCellML normal(p_solver, p_stimulus);
        TS_ASSERT_EQUALS(normal.GetVoltageIndex(), 0u);

        // Optimised model
        Cellluo_rudy_1991FromCellMLOpt opt(p_solver, p_stimulus);
        TS_ASSERT_EQUALS(opt.GetVoltageIndex(), 0u);
        
        // Test the archiving code too
        OutputFileHandler handler("archive", false);
        ArchiveLocationInfo::SetArchiveDirectory(handler.GetOutputDirectoryFullPath());
        std::string archive_filename = ArchiveLocationInfo::GetProcessUniqueFilePath("lr91-pycml.arch");

        // Save both cells at initial state
        {
            AbstractCardiacCell* const p_normal_cell = &normal;
            AbstractCardiacCell* const p_opt_cell = &opt;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_normal_cell;
            output_arch << p_opt_cell;
        }

        //
        // Solve and write to file
        //
        
        // Normal
        ck_start = clock();
        RunOdeSolverWithIonicModel(&normal,
                                   end_time,
                                   "Lr91FromPyCml");
        ck_end = clock();
        double normal_time = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        std::cout << "\n\tNormal: " << normal_time << std::endl;

        CheckCellModelResults("Lr91FromPyCml", "Lr91DelayedStim");

        RunOdeSolverWithIonicModel(&normal,
                                   60.0,
                                   "Lr91GetIIonic", 1000, false);
        TS_ASSERT_DELTA( normal.GetIIonic(), 1.9411, 1e-3);

        // Optimised
        ck_start = clock();
        RunOdeSolverWithIonicModel(&opt,
                                   end_time,
                                   "Lr91FromPyCmlOpt");
        ck_end = clock();
        double opt_time = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        std::cout << "\n\tOptimised: " << opt_time << std::endl;

        CompareCellModelResults("Lr91DelayedStim", "Lr91FromPyCmlOpt", 1e-4, true);

        RunOdeSolverWithIonicModel(&opt,
                                   60.0,
                                   "Lr91GetIIonicOpt", 1000, false);
        TS_ASSERT_DELTA( opt.GetIIonic(), 1.9411, 1e-3);
        
        // Load and check simulation results still match
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractCardiacCell* p_normal_cell;
            AbstractCardiacCell* p_opt_cell;
            input_arch >> p_normal_cell;
            input_arch >> p_opt_cell;

            TS_ASSERT_EQUALS( p_normal_cell->GetNumberOfStateVariables(), 8u );
            TS_ASSERT_EQUALS( p_opt_cell->GetNumberOfStateVariables(), 8u );

            RunOdeSolverWithIonicModel(p_normal_cell,
                                       end_time,
                                       "Lr91FromPyCmlAfterArchive");
            CheckCellModelResults("Lr91FromPyCmlAfterArchive", "Lr91DelayedStim");
            
            RunOdeSolverWithIonicModel(p_opt_cell,
                                       end_time,
                                       "Lr91FromPyCmlOptAfterArchive");
            CompareCellModelResults("Lr91DelayedStim", "Lr91FromPyCmlOptAfterArchive", 1e-4, true);

            delete p_normal_cell;
            delete p_opt_cell;
        }
     }
};


#endif //_TESTPYCML_HPP_
