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


#ifndef _TESTIONICMODELSLONG_HPP_
#define _TESTIONICMODELSLONG_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <ctime>

#include "RunAndCheckIonicModels.hpp"

#include "SimpleStimulus.hpp"
#include "RegularStimulus.hpp"

#include "EulerIvpOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "FoxModel2002Modified.hpp"
#include "BackwardEulerFoxModel2002Modified.hpp"
#include "Maleckar2009OdeSystem.hpp"
#include "Mahajan2008OdeSystem.hpp"
#include "TenTusscher2006OdeSystem.hpp"
#include "CellProperties.hpp"

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
        double duration  = 1.0 ;  // ms
        double start = 50.0; // ms
        double period = 500; // ms
        boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(magnitude, duration, period, start));

        double end_time = 1000.0; //One second in milliseconds


        HeartConfig::Instance()->SetOdeTimeStep(0.002); // 0.005 leads to NaNs.

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        FoxModel2002Modified fox_ode_system(p_solver, p_stimulus);

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
        HeartConfig::Instance()->SetOdeTimeStep(0.01);
        BackwardEulerFoxModel2002Modified backward_system(p_stimulus);
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

    /**
     *  Here we test the scale factors methiods for the mahajan model.
     *  The idea is to set the scale factors for the 3 different cell types (epi, mid and endo)
     *  and check that the rsulting APD makes sense if compared to experiemntal results
     */
    void TestScaleFactorsForMahajanModel(void) throw(Exception)
    {
        double end_time=300;
        double time_step=0.01;
        double sampling_time=time_step;

        // Set stimulus
        double magnitude_stimulus = -70.0;   // pA/pF
        double duration_stimulus = 1.0;  // ms
        double start_stimulus = 10.0;   // ms
        double period=1000;//here, this is ms
        boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(magnitude_stimulus,
                                                                          duration_stimulus,
                                                                          period,
                                                                          start_stimulus));

        boost::shared_ptr<EulerIvpOdeSolver> p_forward_solver(new EulerIvpOdeSolver); //define the solver
        Mahajan2008OdeSystem forward_model(p_forward_solver, p_stimulus);
        boost::shared_ptr<BackwardEulerIvpOdeSolver> p_backward_solver(new BackwardEulerIvpOdeSolver(
                                    forward_model.GetNumberOfStateVariables()));

        Mahajan2008OdeSystem epicardial_model(p_backward_solver, p_stimulus);
        Mahajan2008OdeSystem endocardial_model(p_backward_solver, p_stimulus);
        Mahajan2008OdeSystem midmyocardial_model(p_backward_solver, p_stimulus);

        epicardial_model.SetScaleFactorGks(1.0);
        epicardial_model.SetScaleFactorIto(1.0);
        epicardial_model.SetScaleFactorGkr(1.0);

        midmyocardial_model.SetScaleFactorGks(0.09);
        midmyocardial_model.SetScaleFactorIto(1.0);
        midmyocardial_model.SetScaleFactorGkr(1.0);

        endocardial_model.SetScaleFactorGks(0.86);
        endocardial_model.SetScaleFactorIto(0.2);
        endocardial_model.SetScaleFactorGkr(1.0);

        std::vector<double> state_variables_epi = epicardial_model.GetInitialConditions();
        std::vector<double> state_variables_endo = endocardial_model.GetInitialConditions();
        std::vector<double> state_variables_mid = midmyocardial_model.GetInitialConditions();

        const std::string mahajan_epi_file = "Mahajan_epi";
        const std::string mahajan_mid_file = "Mahajan_mid";
        const std::string mahajan_endo_file = "Mahajan_endo";

        // Solve and write to file

        OdeSolution epi_solution;
        epi_solution = p_backward_solver->Solve(&epicardial_model, state_variables_epi, 0, end_time, time_step, sampling_time);

        epi_solution.WriteToFile("TestIonicModels",
                                 mahajan_epi_file,
                                 &epicardial_model,
                                 "ms",//time units
                                 100,//steps per row
                                 false);/*true cleans the directory*/

        OdeSolution mid_solution;
        mid_solution = p_backward_solver->Solve(&midmyocardial_model, state_variables_mid, 0, end_time, time_step, sampling_time);

        mid_solution.WriteToFile("TestIonicModels",
                                 mahajan_mid_file,
                                 &epicardial_model,
                                 "ms",//time units
                                 100,//steps per row
                                 false);/*true cleans the directory*/

        OdeSolution endo_solution;
        endo_solution = p_backward_solver->Solve(&endocardial_model, state_variables_endo, 0, end_time, time_step, sampling_time);

        endo_solution.WriteToFile("TestIonicModels",
                                  mahajan_endo_file,
                                  &midmyocardial_model,
                                  "ms",//time units
                                  100,//steps per row
                                  false);/*true cleans the directory*/


        ColumnDataReader data_reader_epi("TestIonicModels", mahajan_epi_file);
        ColumnDataReader data_reader_mid("TestIonicModels", mahajan_mid_file);
        ColumnDataReader data_reader_endo("TestIonicModels", mahajan_endo_file);

        std::vector<double> times = data_reader_epi.GetValues("Time");
        std::vector<double> v_endo = data_reader_endo.GetValues("V");
        std::vector<double> v_epi = data_reader_epi.GetValues("V");
        std::vector<double> v_mid = data_reader_mid.GetValues("V");

        CellProperties  cell_properties_endo(v_endo, times);
        CellProperties  cell_properties_epi(v_epi, times);
        CellProperties  cell_properties_mid(v_mid, times);

        double epi_APD = cell_properties_epi.GetLastActionPotentialDuration(90);
        double endo_APD = cell_properties_endo.GetLastActionPotentialDuration(90);
        double mid_APD = cell_properties_mid.GetLastActionPotentialDuration(90);

        //check that percentage increase from epi to mid and endo (roughly) matches results
        // from McIntosh et al. Card Res, 45:397-409. 200 (Figure 1 and 2)
        TS_ASSERT_DELTA((mid_APD-epi_APD)*100/epi_APD, 36.2, 2);
        TS_ASSERT_DELTA((endo_APD-epi_APD)*100/epi_APD, 8, 2);

    }
    
    /**
     * Here we test that the scale factors for the TT model do what they are expected to
     * We check that they modify APD in a way that is expected.
     */
    void TestScaleFactorsForTT06(void)
    {  
        double simulation_end=500;/*end time, in milliseconds for this model*/

        // Set the stimulus, the following values are appropriate for single cell simulations of this model.
        double magnitude = -38.0;   // pA/pF
        double duration = 1.0;  // ms
        double start = 100;   // ms
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude,
                                                                        duration,
                                                                        start));
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver); //define the solver
        HeartConfig::Instance()->SetOdeTimeStep(0.001);// with Forward Euler, this must be as small as 0.001.

        
        const std::string control_file = "TT_epi";
        const std::string mid_file = "TT_mid";
        const std::string endo_file = "TT_endo";   
        const std::string LQT_file = "TT_LQT";                                                        

        TenTusscher2006OdeSystem TT_model_epi(p_solver, p_stimulus);
        TT_model_epi.SetScaleFactorGks(1.0);
        TT_model_epi.SetScaleFactorIto(1.0);
        TT_model_epi.SetScaleFactorGkr(1.0);
        //run the model
        RunOdeSolverWithIonicModel(&TT_model_epi,
                                   simulation_end,
                                   control_file,
                                   100,
                                   false);
                                   
        TenTusscher2006OdeSystem TT_model_mid(p_solver, p_stimulus);                     
        TT_model_mid.SetScaleFactorGks(0.25);
        TT_model_mid.SetScaleFactorIto(1.0);
        TT_model_mid.SetScaleFactorGkr(1.0);

        RunOdeSolverWithIonicModel(&TT_model_mid,
                                   simulation_end,
                                   mid_file,
                                   100,
                                   false);
                                   
        TenTusscher2006OdeSystem TT_model_endo(p_solver, p_stimulus); 
                                  
        TT_model_endo.SetScaleFactorGks(0.66);
        TT_model_endo.SetScaleFactorIto(0.165);
        TT_model_endo.SetScaleFactorGkr(1.0);

        RunOdeSolverWithIonicModel(&TT_model_endo,
                                   simulation_end,
                                   endo_file,
                                   100,
                                   false);
        
        TenTusscher2006OdeSystem TT_model_LQT(p_solver, p_stimulus);                            
        TT_model_LQT.SetScaleFactorGks(1.0);
        TT_model_LQT.SetScaleFactorIto(1.0);
        TT_model_LQT.SetScaleFactorGkr(0.0);

        RunOdeSolverWithIonicModel(&TT_model_LQT,
                                   simulation_end,
                                   LQT_file,
                                   100,
                                   false);

        ColumnDataReader data_reader1("TestIonicModels", control_file);
        std::vector<double> voltages1 = data_reader1.GetValues("V");
        ColumnDataReader data_reader2("TestIonicModels", mid_file);
        std::vector<double> voltages2 = data_reader2.GetValues("V");
        ColumnDataReader data_reader3("TestIonicModels", endo_file);
        std::vector<double> voltages3 = data_reader3.GetValues("V");
        ColumnDataReader data_reader4("TestIonicModels", LQT_file);
        std::vector<double> voltages4 = data_reader4.GetValues("V");
        
        TS_ASSERT_EQUALS(voltages1.size(), voltages2.size());
        TS_ASSERT_EQUALS(voltages2.size(), voltages3.size());
        TS_ASSERT_EQUALS(voltages3.size(), voltages4.size());
        
        //create the times vector
        std::vector<double> times;
        double k =0;
        for (unsigned i=0; i<voltages2.size(); i++)
        {
          times.push_back(k);
          k=k+0.1;
        }     
        
        CellProperties  cell_properties_control(voltages1, times);
        CellProperties  cell_properties_mid(voltages2, times);
        CellProperties  cell_properties_endo(voltages3, times);
        CellProperties  cell_properties_LQT(voltages4, times);
        
        double control_APD = cell_properties_control.GetLastActionPotentialDuration(90);
        double mid_APD = cell_properties_mid.GetLastActionPotentialDuration(90);
        double endo_APD = cell_properties_endo.GetLastActionPotentialDuration(90);
        double LQT_APD = cell_properties_LQT.GetLastActionPotentialDuration(90);

        TS_ASSERT_DELTA(control_APD, 299.82, 0.1);
        TS_ASSERT_DELTA(mid_APD, 391.52, 0.1);
        TS_ASSERT_DELTA(endo_APD, 328.52, 0.1);
        TS_ASSERT_DELTA(LQT_APD , 347.125, 0.1);
     }
         
    void TestScaleFactorsMaleckar(void) throw (Exception)
    {
        double end_time =500;
        boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(-280,
                                                                          6,
                                                                          1000,
                                                                          4.0));

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver); //define the solver
        HeartConfig::Instance()->SetOdeTimeStep(0.001);
        Maleckar2009OdeSystem atrial_ode_system(p_solver, p_stimulus);
        
        const std::string control_file = "control";
        const std::string first_set_file = "first_scale_factor_set";
        const std::string second_set_file = "second_scale_factor_set";
        const std::string AZD_file = "AZD_scale_factor_set";
        
        OdeSolution control_solution;
        OdeSolution first_scale_factor_set_solution;
        OdeSolution second_scale_factor_set_solution;
        OdeSolution AZD_scale_factor_set_solution;
        
        double time_step=0.001;
        double sampling_time=0.001;       
        std::vector<double> state_variables= atrial_ode_system.GetInitialConditions();
        
        //default values
        atrial_ode_system.SetScaleFactorGks(1.0);
        atrial_ode_system.SetScaleFactorIto(1.0);
        atrial_ode_system.SetScaleFactorGkr(1.0);
        atrial_ode_system.SetScaleFactorGna(1.0);
        atrial_ode_system.SetScaleFactorAch(1e-24);
        atrial_ode_system.SetScaleFactorGNaK(1.0);
        atrial_ode_system.SetScaleFactorGNaCa(1.0);
        atrial_ode_system.SetScaleFactorGKur(1.0);
        atrial_ode_system.SetScaleFactorGK1(1.0);
        atrial_ode_system.SetScaleFactorGCaL(1.0);
        atrial_ode_system.SetScaleFactorAZD(0.0);
            
        control_solution = p_solver->Solve(&atrial_ode_system, state_variables, 0, end_time, time_step, sampling_time);     
        control_solution.WriteToFile("TestIonicModels",
                              control_file,
                              &atrial_ode_system,
                              "ms",//time units
                              100,//steps per row
                              false);/*true cleans the directory*/
                              
        //now apply the first scale factor set, decreases outward currents                      
        atrial_ode_system.SetScaleFactorGks(0.8);
        atrial_ode_system.SetScaleFactorIto(1.0);
        atrial_ode_system.SetScaleFactorGkr(0.9);
        atrial_ode_system.SetScaleFactorGna(1.0);
        atrial_ode_system.SetScaleFactorAch(1e-24);
        atrial_ode_system.SetScaleFactorGNaK(1.0);
        atrial_ode_system.SetScaleFactorGNaCa(1.0);
        atrial_ode_system.SetScaleFactorGKur(0.7);
        atrial_ode_system.SetScaleFactorGK1(0.6);
        atrial_ode_system.SetScaleFactorGCaL(1.0);
        atrial_ode_system.SetScaleFactorAZD(0.0);
        
        state_variables= atrial_ode_system.GetInitialConditions();
        first_scale_factor_set_solution = p_solver->Solve(&atrial_ode_system, state_variables, 0, end_time, time_step, sampling_time);            
        first_scale_factor_set_solution.WriteToFile("TestIonicModels",
                              first_set_file,
                              &atrial_ode_system,
                              "ms",//time units
                              100,//steps per row
                              false);/*true cleans the directory*/
         
        //now apply the secondscale factor set, this one increases inward currents                  
        atrial_ode_system.SetScaleFactorGks(1.0);
        atrial_ode_system.SetScaleFactorIto(1.0);
        atrial_ode_system.SetScaleFactorGkr(1.0);
        atrial_ode_system.SetScaleFactorGna(1.5);
        atrial_ode_system.SetScaleFactorAch(1e-24);
        atrial_ode_system.SetScaleFactorGNaK(1.0);
        atrial_ode_system.SetScaleFactorGNaCa(2.0);
        atrial_ode_system.SetScaleFactorGKur(1.0);
        atrial_ode_system.SetScaleFactorGK1(1.0);
        atrial_ode_system.SetScaleFactorGCaL(1.6);
        atrial_ode_system.SetScaleFactorAZD(0.0);
        
        state_variables= atrial_ode_system.GetInitialConditions();
        second_scale_factor_set_solution = p_solver->Solve(&atrial_ode_system, state_variables, 0, end_time, time_step, sampling_time);            
        second_scale_factor_set_solution.WriteToFile("TestIonicModels",
                              second_set_file,
                              &atrial_ode_system,
                              "ms",//time units
                              100,//steps per row
                              false);/*true cleans the directory*/
        
        //check the AZD scale factor (vs control)
        atrial_ode_system.SetScaleFactorGks(1.0);
        atrial_ode_system.SetScaleFactorIto(1.0);
        atrial_ode_system.SetScaleFactorGkr(1.0);
        atrial_ode_system.SetScaleFactorGna(1.0);
        atrial_ode_system.SetScaleFactorAch(1e-24);
        atrial_ode_system.SetScaleFactorGNaK(1.0);
        atrial_ode_system.SetScaleFactorGNaCa(1.0);
        atrial_ode_system.SetScaleFactorGKur(1.0);
        atrial_ode_system.SetScaleFactorGK1(1.0);
        atrial_ode_system.SetScaleFactorGCaL(1.0);
        atrial_ode_system.SetScaleFactorAZD(5.0);
            
        AZD_scale_factor_set_solution = p_solver->Solve(&atrial_ode_system, state_variables, 0, end_time, time_step, sampling_time);     
        AZD_scale_factor_set_solution.WriteToFile("TestIonicModels",
                              AZD_file,
                              &atrial_ode_system,
                              "ms",//time units
                              100,//steps per row
                              false);/*true cleans the directory*/
                              
                                                       
        ColumnDataReader data_reader1("TestIonicModels", control_file);
        std::vector<double> voltages1 = data_reader1.GetValues("V");
        ColumnDataReader data_reader2("TestIonicModels", first_set_file);
        std::vector<double> voltages2 = data_reader2.GetValues("V");
        ColumnDataReader data_reader3("TestIonicModels", second_set_file);
        std::vector<double> voltages3 = data_reader3.GetValues("V");
        ColumnDataReader data_reader4("TestIonicModels", AZD_file);
        std::vector<double> voltages4 = data_reader4.GetValues("V");
        
        TS_ASSERT_EQUALS(voltages1.size(), voltages2.size());
        TS_ASSERT_EQUALS(voltages2.size(), voltages3.size());
        TS_ASSERT_EQUALS(voltages3.size(), voltages4.size());
        
        //create the times vector
        std::vector<double> times;
        double k =0;
        for (unsigned i=0; i<voltages2.size(); i++)
        {
          times.push_back(k);
          k=k+0.1;
        }     
        
        CellProperties  cell_properties_control(voltages1, times);
        CellProperties  cell_properties_first(voltages2, times);
        CellProperties  cell_properties_second(voltages3, times);
        CellProperties  cell_properties_AZD(voltages4, times);
        
        double control_APD = cell_properties_control.GetLastActionPotentialDuration(90);
        double first_APD = cell_properties_first.GetLastActionPotentialDuration(90);
        double second_APD = cell_properties_second.GetLastActionPotentialDuration(90);
        double AZD_APD = cell_properties_AZD.GetLastActionPotentialDuration(90);
        
        //test that the aps are actually longer than control (all interventions were meant to have that effect except the last one)
        TS_ASSERT_LESS_THAN(control_APD, first_APD);
        TS_ASSERT_LESS_THAN(control_APD, second_APD);
        TS_ASSERT_LESS_THAN(AZD_APD, control_APD);
        
        //leave some hardcoded value for testing
        TS_ASSERT_DELTA(control_APD, 200.36, 0.1);
        TS_ASSERT_DELTA(first_APD, 384.77, 0.1);
        TS_ASSERT_DELTA(second_APD, 308.3515, 0.1);
        TS_ASSERT_DELTA(AZD_APD , 183.67, 0.1);
    }
};


#endif //_TESTIONICMODELSLONG_HPP_
