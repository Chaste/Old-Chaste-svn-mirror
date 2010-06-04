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
#include "VectorHelperFunctions.hpp"

#include "luo_rudy_1991.hpp"
#include "luo_rudy_1991Opt.hpp"
#include "luo_rudy_1991BackwardEuler.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"

#ifdef CHASTE_CVODE
#include "luo_rudy_1991Cvode.hpp"
#include "luo_rudy_1991CvodeOpt.hpp"
#endif // CHASTE_CVODE

class TestPyCml : public CxxTest::TestSuite
{
    template<typename VECTOR_TYPE>
    void CheckDerivedQuantities(AbstractParameterisedSystem<VECTOR_TYPE>& rCell,
                                const VECTOR_TYPE& rStateVec)
    {
        TS_ASSERT_EQUALS(rCell.GetNumberOfDerivedQuantities(), 2u);
        TS_ASSERT_EQUALS(rCell.GetDerivedQuantityIndex("FonRT"), 0u);
        TS_ASSERT_EQUALS(rCell.GetDerivedQuantityIndex("potassium_currents"), 1u);
        TS_ASSERT_EQUALS(rCell.GetDerivedQuantityUnits(0u), "per_millivolt");
        TS_ASSERT_EQUALS(rCell.GetDerivedQuantityUnits(1u), "microA_per_cm2");
        VECTOR_TYPE derived = rCell.ComputeDerivedQuantitiesFromCurrentState(0.0);
        const double FonRT = 0.037435728309031795;
        const double i_K_total = 1.0007;
        TS_ASSERT_EQUALS(GetVectorSize(derived), 2u);
        TS_ASSERT_DELTA(GetVectorComponent(derived, 0), FonRT, 1e-12);
        TS_ASSERT_DELTA(GetVectorComponent(derived, 1), i_K_total, 1e-4);
        DeleteVector(derived);
        derived = rCell.ComputeDerivedQuantities(0.0, rStateVec);
        TS_ASSERT_EQUALS(GetVectorSize(derived), 2u);
        TS_ASSERT_DELTA(GetVectorComponent(derived, 0), FonRT, 1e-12);
        TS_ASSERT_DELTA(GetVectorComponent(derived, 1), i_K_total, 1e-4);
        DeleteVector(derived);
    }
    
    template<typename VECTOR_TYPE>
    void CheckParameter(AbstractParameterisedSystem<VECTOR_TYPE>& rCell)
    {
        TS_ASSERT_EQUALS(rCell.GetNumberOfParameters(), 1u);
        TS_ASSERT_EQUALS(rCell.GetParameterIndex("fast_sodium_current_conductance"), 0u);
        TS_ASSERT_EQUALS(rCell.GetParameterUnits(0u), "milliS_per_cm2");
        TS_ASSERT_EQUALS(rCell.GetParameter(0u), 23.0);
        rCell.SetParameter(0u, 0.1);
        TS_ASSERT_EQUALS(rCell.GetParameter(0u), 0.1);
        rCell.SetParameter(0u, 23.0);
    }
    
public:
    /** For comparison with the test below; copied from TestIonicModels.hpp */
    void TestOdeSolverForLR91WithDelayedSimpleStimulus(void)
    {
        clock_t ck_start, ck_end;

        // Set stimulus
        double magnitude = -25.5;
        double duration  = 2.0 ;  // ms
        double when = 50.0; // ms

        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude, duration, when));
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        double end_time = 1000.0; //One second in milliseconds

        LuoRudyIModel1991OdeSystem lr91_ode_system(p_solver, p_stimulus);
        TS_ASSERT_EQUALS(lr91_ode_system.GetVoltageIndex(), 4u); // For coverage

        // Solve and write to file
        ck_start = clock();
        RunOdeSolverWithIonicModel(&lr91_ode_system,
                                   end_time,
                                   "Lr91DelayedStim",
                                   100,
                                   false);
        ck_end = clock();
        double forward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        std::cout << "\n\tForward: " << forward << std::endl;
        CheckCellModelResults("Lr91DelayedStim");
    }

    /**
     * This test is designed to quickly check that PyCml-generated code matches the Chaste interfaces,
     * and gives expected results.
     */
    void TestPyCmlCodeGeneration() throw(Exception)
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
        double i_ionic_end_time = 60.0; // ms
        double i_ionic = 1.9411; // test value

        // Normal model
        Cellluo_rudy_1991FromCellML normal(p_solver, p_stimulus);
        TS_ASSERT_EQUALS(normal.GetVoltageIndex(), 0u);

        // Optimised model
        Cellluo_rudy_1991FromCellMLOpt opt(p_solver, p_stimulus);
        TS_ASSERT_EQUALS(opt.GetVoltageIndex(), 0u);
        
        // Backward Euler optimised model
        Cellluo_rudy_1991FromCellMLBackwardEuler be(p_solver, p_stimulus);
        TS_ASSERT_EQUALS(be.GetVoltageIndex(), 0u);

        // Check that the tables exist!
        double v = opt.GetVoltage();
        opt.SetVoltage(-100000);
        TS_ASSERT_THROWS_CONTAINS(opt.GetIIonic(), "V outside lookup table range");
        opt.SetVoltage(v);
        
        be.SetVoltage(-100000);
        TS_ASSERT_THROWS_CONTAINS(be.GetIIonic(), "V outside lookup table range");
        be.SetVoltage(v);
        
        // Single parameter
        CheckParameter(normal);
        CheckParameter(opt);
        CheckParameter(be);
        
        // Derived variables
        CheckDerivedQuantities(normal, normal.GetInitialConditions());
        CheckDerivedQuantities(opt, opt.GetInitialConditions());
        CheckDerivedQuantities(be, be.GetInitialConditions());

#ifdef CHASTE_CVODE
        // CVODE version
        Cellluo_rudy_1991FromCellMLCvode cvode_cell(p_solver, p_stimulus);
        TS_ASSERT_EQUALS(cvode_cell.GetVoltageIndex(), 0u);
        // Optimised CVODE version
        Cellluo_rudy_1991FromCellMLCvodeOpt cvode_opt(p_solver, p_stimulus);
        TS_ASSERT_EQUALS(cvode_opt.GetVoltageIndex(), 0u);
        // Check that the tables exist!
        v = opt.GetVoltage();
        cvode_opt.SetVoltage(-100000);
        TS_ASSERT_THROWS_CONTAINS(cvode_opt.GetIIonic(), "V outside lookup table range");
        cvode_opt.SetVoltage(v);
        
        // Single parameter
        CheckParameter(cvode_cell);
        CheckParameter(cvode_opt);
        
        // Derived variables
        N_Vector inits = cvode_cell.GetInitialConditions();
        CheckDerivedQuantities(cvode_cell, inits);
        CheckDerivedQuantities(cvode_opt, inits);
        DeleteVector(inits);
#endif // CHASTE_CVODE

        // Test the archiving code too
        OutputFileHandler handler("archive", false);
        handler.SetArchiveDirectory();
        std::string archive_filename = ArchiveLocationInfo::GetProcessUniqueFilePath("lr91-pycml.arch");

        // Save all (non-CVODE) cells at initial state
        {
            AbstractCardiacCell* const p_normal_cell = &normal;
            AbstractCardiacCell* const p_opt_cell = &opt;
            AbstractCardiacCell* const p_be_cell = &be;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_normal_cell;
            output_arch << p_opt_cell;
            output_arch << p_be_cell;
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
                                   i_ionic_end_time,
                                   "Lr91GetIIonic", 1000, false);
        TS_ASSERT_DELTA( normal.GetIIonic(), i_ionic, 1e-3);
        
        // With zero g_Na
        normal.SetParameter(0u, 0.0);
        normal.SetStateVariables(normal.GetInitialConditions());
        RunOdeSolverWithIonicModel(&normal,
                                   end_time,
                                   "Lr91FromPyCmlZeroGna");
        CheckCellModelResults("Lr91FromPyCmlZeroGna");

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
                                   i_ionic_end_time,
                                   "Lr91GetIIonicOpt", 1000, false);
        TS_ASSERT_DELTA( opt.GetIIonic(), i_ionic, 1e-3);
        
        // No stimulus at end time
        TS_ASSERT_DELTA(opt.Get_membrane__I_stim(), 0.0, 1e-12);
        
        // Backward Euler
        ck_start = clock();
        RunOdeSolverWithIonicModel(&be,
                                   end_time,
                                   "Lr91FromPyCmlBackwardEuler");
        ck_end = clock();
        double be_time = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        std::cout << "\n\tBackward Euler: " << be_time << std::endl;

        CompareCellModelResults("Lr91DelayedStim", "Lr91FromPyCmlBackwardEuler", 1e-2, true);

        RunOdeSolverWithIonicModel(&be,
                                   i_ionic_end_time,
                                   "Lr91GetIIonicBackwardEuler", 1000, false);
        TS_ASSERT_DELTA( be.GetIIonic(), i_ionic, 1e-3);

#ifdef CHASTE_CVODE
        // CVODE
        double max_dt = 1.0; //ms
        ck_start = clock();
        OdeSolution cvode_solution = cvode_cell.Solve(0.0, end_time, max_dt, max_dt);
        cvode_solution.WriteToFile("TestIonicModels","Lr91FromPyCmlCvode","ms",1,false);
        ck_end = clock();
        double cvode_time = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        std::cout << "\n\tCVODE: " << cvode_time << std::endl;
        CompareCellModelResults("Lr91FromPyCml", "Lr91FromPyCmlCvode", 1e-1, true);
        // Coverage
        cvode_cell.SetVoltageDerivativeToZero();
        cvode_cell.SetStateVariables(cvode_cell.GetInitialConditions());
        cvode_cell.Solve(0.0, i_ionic_end_time, max_dt);
        TS_ASSERT_DELTA(cvode_cell.GetIIonic(), 0.0, 1e-1); // Cell should be at rest
        cvode_cell.SetVoltageDerivativeToZero(false);
        // Check GetIIonic
        cvode_cell.SetStateVariables(cvode_cell.GetInitialConditions());
        cvode_cell.Solve(0.0, i_ionic_end_time, max_dt);
        TS_ASSERT_DELTA(cvode_cell.GetIIonic(), i_ionic, 1e-1);

        // CVODE Optimised
        ck_start = clock();
        cvode_solution = cvode_opt.Solve(0.0, end_time, max_dt, max_dt);
        cvode_solution.WriteToFile("TestIonicModels","Lr91FromPyCmlCvodeOpt","ms",1,false);
        ck_end = clock();
        double cvode_opt_time = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        std::cout << "\n\tCVODE Optimised: " << cvode_opt_time << std::endl;
        CompareCellModelResults("Lr91FromPyCmlCvode", "Lr91FromPyCmlCvodeOpt", 1e-1, true);
        // Coverage
        cvode_opt.SetVoltageDerivativeToZero();
        cvode_opt.SetStateVariables(cvode_opt.GetInitialConditions());
        cvode_opt.Solve(0.0, i_ionic_end_time, max_dt);
        TS_ASSERT_DELTA(cvode_opt.GetIIonic(), 0.0, 1e-1); // Cell should be at rest
        cvode_opt.SetVoltageDerivativeToZero(false);
        // Check GetIIonic
        cvode_opt.SetStateVariables(cvode_opt.GetInitialConditions());
        cvode_opt.Solve(0.0, i_ionic_end_time, max_dt);
        TS_ASSERT_DELTA(cvode_opt.GetIIonic(), i_ionic, 1e-1);
        
        // No stimulus at end time
        TS_ASSERT_DELTA(cvode_opt.Get_membrane__I_stim(), 0.0, 1e-12);
#endif // CHASTE_CVODE

        // Load and check simulation results still match
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractCardiacCell* p_normal_cell;
            AbstractCardiacCell* p_opt_cell;
            AbstractCardiacCell* p_be_cell;
            input_arch >> p_normal_cell;
            input_arch >> p_opt_cell;
            input_arch >> p_be_cell;

            TS_ASSERT_EQUALS( p_normal_cell->GetNumberOfStateVariables(), 8u );
            TS_ASSERT_EQUALS( p_opt_cell->GetNumberOfStateVariables(), 8u );
            TS_ASSERT_EQUALS( p_be_cell->GetNumberOfStateVariables(), 8u );

            RunOdeSolverWithIonicModel(p_normal_cell,
                                       end_time,
                                       "Lr91FromPyCmlAfterArchive");
            CheckCellModelResults("Lr91FromPyCmlAfterArchive", "Lr91DelayedStim");

            RunOdeSolverWithIonicModel(p_opt_cell,
                                       end_time,
                                       "Lr91FromPyCmlOptAfterArchive");
            CompareCellModelResults("Lr91DelayedStim", "Lr91FromPyCmlOptAfterArchive", 1e-4, true);

            RunOdeSolverWithIonicModel(p_be_cell,
                                       end_time,
                                       "Lr91FromPyCmlBackwardEulerAfterArchive");
            CompareCellModelResults("Lr91DelayedStim", "Lr91FromPyCmlBackwardEulerAfterArchive", 1e-2, true);

            delete p_normal_cell;
            delete p_opt_cell;
            delete p_be_cell;
        }
     }
};


#endif //_TESTPYCML_HPP_
