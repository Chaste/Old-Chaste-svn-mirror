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
#ifdef CHASTE_CVODE

#ifndef _TESTMODIFYINGCONDUCTANCES_HPP_
#define _TESTMODIFYINGCONDUCTANCES_HPP_

#include <cxxtest/TestSuite.h>
#include <ctime>
#include <iostream>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/shared_ptr.hpp>

#include "AbstractCardiacCell.hpp"
#include "AbstractCvodeCell.hpp"
#include "AbstractCardiacCellWithModifiers.hpp"
#include "Exception.hpp"
#include "OutputFileHandler.hpp"
#include "Modifiers.hpp"

#include "DataStructure.hpp"
#include "CellProperties.hpp"

#include "S1S2Stimulus.hpp"
#include "DynamicRestitutionStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"

#include "Shannon2004Cvode.hpp"
#include "ten_tusscher_model_2006_epiCvode.hpp"
#include "hund_rudy_2004_aCvode.hpp"
#include "mahajan_2008Cvode.hpp"
#include "grandi2010ssCvode.hpp"

class TestModifyingConductances : public CxxTest::TestSuite
{
public:
    /**
     * This test will wipe $CHASTE_TEST_OUTPUT/VaryingConductances/<Model>/
     *
     * Parameters can be defined at the top of this Test
     */
    void TestDrugAffectByVaryingConductances(void) throw (Exception)
    {
        //////////// DEFINE PARAMETERS ///////////////
        CommandLineArguments* p_args = CommandLineArguments::Instance();
        unsigned argc = *(p_args->p_argc); // has the number of arguments.
        std::cout << "# " << argc-1 << " arguments supplied.\n" << std::flush;

        if (argc == 1)
        {
            std::cerr << "TestModifyingConductances::Please input an argument\n"
                         "* --model\n"
                         "    options: 1= Shannon, 2 = TenTusscher, 3 = Mahajan,\n"
                         "             4 = HundRudy, 5 = Grandi.\n"
                         "* --herg_only\n"
                         "    optional: ignore the Na and CaL channel block data.\n";
            return;
        }

        bool just_herg_block = CommandLineArguments::Instance()->OptionExists("--herg_only");
        char* val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--model");
        unsigned model_index = atoi(val);

        // Compile-time arguments
        bool steady_state_AP_experiment = true;
        bool s1s2_protocol_experiment = true;
        bool dynamic_protocol_experiment = true;

        ///////// END DEFINE PARAMETERS ////////////////////////

        /// Cvode cells use a CVODE solver regardless of which standard solver is passed in.
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        //boost::shared_ptr<CvodeAdaptor> p_solver(new CvodeAdaptor);

        // Load drug data
        DataStructure drug_data("./projects/CardiovascRes11/test/drug_data.dat");

        // Pick the model...
        // 1= Shannon, 2=TenTusscher, 3=Priebe, 4=Livshitz, 5 = Mahajan, 6 = HundRudy
        AbstractCardiacCellWithModifiers<AbstractCvodeCell>* p_model;
        boost::shared_ptr<AbstractStimulusFunction> p_stimulus; // This will then be set to any of the above stimuli in the individual protocols.

        switch (model_index)
        {
            case 1u:
                p_model = new CellShannon2004FromCellMLCvode(p_solver, p_stimulus);
                break;
            case 2u:
                p_model = new Cellten_tusscher_model_2006_epiFromCellMLCvode(p_solver, p_stimulus);
                break;
            case 3u:
                p_model = new Cellmahajan_2008FromCellMLCvode(p_solver, p_stimulus);
                break;
            case 4u:
                p_model = new Cellhund_rudy_2004_aFromCellMLCvode(p_solver, p_stimulus);
                break;
            case 5u:
                p_model = new Cellgrandi2010ssFromCellMLCvode(p_solver, p_stimulus);
                break;
            default:
                EXCEPTION("No model matches this index");
        }

        double s_magnitude = -15; // We will attempt to overwrite these with model specific ones below
        double s_duration = 3.0; // We will attempt to overwrite these with model specific ones below
        // Use the default CellML stimulus amplitude and duration, but set start time and period to what we want.

        if (p_model->HasCellMLDefaultStimulus())
        {
            p_model->UseCellMLDefaultStimulus();
            boost::shared_ptr<RegularStimulus> p_reg_stim =
                     boost::static_pointer_cast<RegularStimulus>(p_model->GetStimulusFunction());
            s_magnitude = p_reg_stim->GetMagnitude();
            s_duration = p_reg_stim->GetDuration();
        }

        double s_start = 1.0;           // ms
        double hertz = 1.0;             // s^-1
        double s1_period = 1000/hertz;  // ms
        double printing_time_step = 1.0; // ms
        double maximum_time_step; // ms
        if (printing_time_step > s_duration)
        {
            maximum_time_step = s_duration;
        }
        else
        {
            maximum_time_step = printing_time_step;
        }


        // PARAMETERS FOR THE STEADY STATE EXPERIMENTS
        // This is a regular stimulus for steady-state pacing experiments
        boost::shared_ptr<RegularStimulus> p_regular_stimulus(new RegularStimulus(s_magnitude, s_duration, s1_period, s_start));
        const unsigned num_stims = 1000u;

        // PARAMETERES FOR THE S1S2 PROTOCOL //
        double duration_of_s1 = num_stims*s1_period; // 100 APs
        double s2_pcls[28] = {1000, 950, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400,
                     380, 360, 340, 320, 300, 280, 260, 240, 220, 200, 180, 160, 140, 120, 100};
        std::vector<double> s2_periods(s2_pcls,s2_pcls+28);
        boost::shared_ptr<S1S2Stimulus> p_s1s2_stimulus(new S1S2Stimulus(s_magnitude, s_duration, duration_of_s1, s1_period, s_start, s2_periods));

        // PARAMETERS FOR THE DYNAMIC RESTITUTION PROTOCOL //
        double dynamic_pcls[60] = {1000, 900, 800, 700, 600, 550, 500, 480, 460, 440, 420, 400,
                380, 360, 350, 340, 330, 320, 310, 300, 295, 290, 285, 280, 275, 270, 265, 260,
                255, 250, 245, 240, 235, 230, 225, 220, 215, 210, 205, 200, 195, 190, 185, 180,
                175, 170, 165, 160, 155, 150, 145, 140, 135, 130, 125, 120, 115, 110, 105, 100};
        std::vector<double> pacing_cycle_lengths(dynamic_pcls,dynamic_pcls+60);
        const unsigned num_pulses = 100u;
        const unsigned num_pulses_to_analyse = 8;
        boost::shared_ptr<DynamicRestitutionStimulus> p_dynamic_stimulus(new DynamicRestitutionStimulus(s_magnitude, s_duration, s_start, pacing_cycle_lengths, num_pulses));

        boost::shared_ptr<const AbstractOdeSystemInformation> p_ode_info = p_model->GetSystemInformation();
        std::string model_name = p_model->GetSystemName();
        unsigned calcium_index = p_model->GetAnyVariableIndex("cytosolic_calcium_concentration");
        unsigned voltage_index = p_ode_info->GetStateVariableIndex("membrane_voltage");

        // Set up foldernames for each model and protocol set.
        std::string foldername = "VaryingConductances/" + model_name;
        if (just_herg_block)
        {
            foldername = foldername + "HergOnly";
        }
        std::stringstream steady_state_foldername;
        std::stringstream s1s2_foldername;
        std::stringstream dynamic_foldername;

        steady_state_foldername << foldername << "/" << hertz << "Hz_SteadyState";
        s1s2_foldername << foldername << "/S1S2";
        dynamic_foldername << foldername << "/Dynamic";

        // Make and clean the above directories.
        OutputFileHandler steady_handler(steady_state_foldername.str());
        OutputFileHandler s1s2_handler(s1s2_foldername.str());
        OutputFileHandler dynamic_handler(dynamic_foldername.str());

        // Open files and write headers
        out_stream steady_voltage_results_file =  steady_handler.OpenOutputFile("voltage_results.dat");
        out_stream steady_calcium_results_file =  steady_handler.OpenOutputFile("calcium_results.dat");
        out_stream s1s2_results_file = s1s2_handler.OpenOutputFile("voltage_results.dat");
        out_stream dynamic_results_file = dynamic_handler.OpenOutputFile("voltage_results.dat");
        *steady_voltage_results_file << "DrugIndex\tConcentration(nM)\tUpstrokeVelocity(mV/ms)\tPeakVm(mV)\tAPD50(ms)\tAPD90(ms)\n";
        *steady_calcium_results_file << "DrugIndex\tConcentration(nM)\tPeakChangeInCa(mM/ms)\tPeakCa(mM)\tCaD50(ms)\tCaD90(ms)\n";
        *s1s2_results_file << "DrugIndex\tConcentration(nM)\tS2_PCL(ms)\tAPD50_1(ms)\tAPD90_1(ms)\tAPD50_2(ms)\tAPD90_2(ms)\n";
        *dynamic_results_file << "DrugIndex\tConcentration(nM)\tPCL(ms)";
        for (unsigned i=0; i< num_pulses_to_analyse; i++)
        {
            *dynamic_results_file << "\tAPD90_" << i+1 << "(ms)";
        }
        *dynamic_results_file << "\n";

        /**
         * START LOOP OVER EACH DRUG
         *
         */
        std::cout << "Running simulations..." << std::endl;
        for(unsigned drug_index = 0; drug_index < drug_data.GetNumDrugs(); drug_index++)
        {
            // Loop over each drug
            std::cout << "DRUG = " << drug_data.GetDrugName(drug_index) << "\n" << std::flush;

            // Following code prevents us running analyses if we are missing IC50 values.
            bool three_ic50_values_present = true;
            bool herg_ic50_value_present = true;
            for (unsigned i=0; i<3 ; i++)
            {
                if (drug_data.GetIC50Value(drug_index,i)<0)
                {
                    three_ic50_values_present = false;
                }
            }
            if (drug_data.GetIC50Value(drug_index,2)<0)
            {
                herg_ic50_value_present = false;
            }

            // If we don't have three IC50 values and we are not supposed to be doing herg block
            if ( (three_ic50_values_present==false && !just_herg_block)
                  || herg_ic50_value_present==false )
            {
                continue; // Jump to next drug.
            }

            /**
             * DECIDE WHAT CONCENTRATIONS TO TEST AT
             */
            std::vector<double> drug_conc;
            {
                // This range of drug concentrations corresponds to those shown in Figure 1 of
                // Redfern-et-al (2003) doi:10.1016/S0008-6363(02)00846-5
                drug_conc.push_back(0);
                for (unsigned i=0 ; i<7 ;i++)
                {
                    drug_conc.push_back(pow(10.0,(double)(i)));
                }

                // Add the Clinical dose range for this drug too if that information exists...:
                // again mostly from Fig 1, Redfern-et-al (2003) doi:10.1016/S0008-6363(02)00846-5
                try
                {
                    double low = drug_data.GetClinicalDoseRange(drug_index,0u);
                    double high = drug_data.GetClinicalDoseRange(drug_index,1u);
                    drug_conc.push_back(low);               // low dose
                    drug_conc.push_back((low+high)/2.0);    // mean dose
                    drug_conc.push_back(high);              // high dose
                    drug_conc.push_back(10*high);           // overdose
                }
                catch (Exception &e)
                {   // Check it fell over for the right reason and carry on.
                    TS_ASSERT_EQUALS(e.GetShortMessage(),"No data available on clinical dose for this drug.");
                }
            }

            // Reset to initial conditions (avoid starting from a drugged up state)
            p_model->SetStateVariables(p_model->GetInitialConditions());

            double voltage_AP_threshold; // We will set this uniquely for each model at control dose.

            /**
             * START LOOP OVER EACH CONCENTRATION TO TEST WITH
             */
            for (unsigned conc_index=0; conc_index<drug_conc.size(); conc_index++)
            {
                std::cout << "Drug Conc = " << drug_conc[conc_index] << "\n" << std::flush;

                // Here we calculate the proportion of the different channels which are still active
                // (at this concentration of this drug)
                double gNa_factor = DataStructure::CalculateConductanceFactor(drug_conc[conc_index],drug_data.GetIC50Value(drug_index,0));
                double gCaL_factor = DataStructure::CalculateConductanceFactor(drug_conc[conc_index],drug_data.GetIC50Value(drug_index,1));
                double gKr_factor = DataStructure::CalculateConductanceFactor(drug_conc[conc_index],drug_data.GetIC50Value(drug_index,2));

                if (just_herg_block)
                {
                    gCaL_factor = 1.0;
                    gNa_factor = 1.0;
                }

                std::cout << "gNa factor = " << gNa_factor << "\n" << std::flush;
                std::cout << "gCaL factor = " << gCaL_factor << "\n" << std::flush;
                std::cout << "gKr factor = " << gKr_factor << "\n" << std::flush;

                // 'Factor Modifier' classes allow us to intervene in the right hand side of the ODE,
                // multiplying the usual values of the labelled parameter by a scaling factor.
                boost::shared_ptr<AbstractModifier> a = boost::shared_ptr<AbstractModifier>(new FactorModifier(gNa_factor));
                boost::shared_ptr<AbstractModifier> b = boost::shared_ptr<AbstractModifier>(new FactorModifier(gCaL_factor));
                boost::shared_ptr<AbstractModifier> c = boost::shared_ptr<AbstractModifier>(new FactorModifier(gKr_factor));

                // The following names are fixed and correspond to metadata tags.
                p_model->SetModifier("membrane_fast_sodium_current_conductance",a);
                p_model->SetModifier("membrane_L_type_calcium_current_conductance",b);
                p_model->SetModifier("membrane_rapid_delayed_rectifier_potassium_current_conductance",c);

                /**
                 * STEADY STATE PACING EXPERIMENT
                 */
                if (steady_state_AP_experiment)
                {
                    p_model->SetStimulusFunction(p_regular_stimulus); // Assign the regular stimulus to the cell's stimulus
                    // Stimulate many times to establish the steady state response.
                    for (unsigned j=0 ; j<num_stims; j++)
                    {
                        // Set the max time step to be the printing time step (1ms).
                        double step_to_use = maximum_time_step;
                        if (j==num_stims-1)
                        {
                            step_to_use /= 100; // Get a lot more detail this step.
                        }
                        OdeSolution solution = p_model->Solve(s1_period*(double)j, s1_period*(double)(j+1), step_to_use, step_to_use);
                        if (j==num_stims-1) // If this is the last stimulus write out the results
                        {
                            // Get voltage properties
                            std::vector<double> voltages = solution.GetVariableAtIndex(voltage_index); // Voltage should always be zero

                            if (conc_index==0) // If this is the control case set a sensible threshold for APs
                            {
                                double max_v = -DBL_MAX;
                                double min_v = DBL_MAX;
                                for (unsigned i=0 ; i<voltages.size(); i++)
                                {
                                    if (voltages[i] > max_v)
                                    {
                                        max_v = voltages[i];
                                    }
                                    if (voltages[i] < min_v)
                                    {
                                        min_v = voltages[i];
                                    }
                                }
                                voltage_AP_threshold = min_v + (max_v - min_v)*0.1; // Say threshold for an AP is 10% of way up upstroke.
                            }
                            CellProperties voltage_properties(voltages, solution.rGetTimes(),voltage_AP_threshold);// Threshold for AP is above -70ish.

                            double apd90, apd50, upstroke, peak;
                            try
                            {
                                apd90 = voltage_properties.GetLastActionPotentialDuration(90);
                                apd50 = voltage_properties.GetLastActionPotentialDuration(50);
                                upstroke = voltage_properties.GetLastCompleteMaxUpstrokeVelocity();
                                peak = voltage_properties.GetLastCompletePeakPotential();
                            }
                            catch (Exception& e)
                            {
                                if (e.GetShortMessage()=="AP did not occur, never exceeded threshold voltage." ||
                                        e.GetShortMessage()=="No full action potential was recorded")
                                {
                                    std::cout << "At Steady pacing of " << hertz << "Hz we did not record any APs\n" << std::flush;
                                    apd90 = -1;
                                    apd50 = -1;
                                    upstroke = -1;
                                    peak = -1; // Default error/missing data code for matlab analysis.
                                }
                                else
                                {
                                    throw e;
                                }
                            }

                            std::cout << "1Hz Upstroke velocity = " << upstroke << ", Peak mV = " << peak << ", APD50 = " << apd50 << ", APD90 = " << apd90 << "\n" << std::flush;
                            *steady_voltage_results_file << drug_index << "\t" << drug_conc[conc_index] << "\t" << upstroke << "\t" << peak << "\t" << apd50  << "\t" << apd90 << "\n";

                            // Get Calcium properties
                            std::vector<double> calcium = solution.GetVariableAtIndex(calcium_index);
                            // Calculate threshold based on max - min of calcium transient
                            double max_ca = -DBL_MAX;
                            double min_ca = DBL_MAX;
                            for (unsigned i=0 ; i<calcium.size(); i++)
                            {
                                if (calcium[i] > max_ca)
                                {
                                    max_ca = calcium[i];
                                }
                                if (calcium[i] < min_ca)
                                {
                                    min_ca = calcium[i];
                                }
                            }
                            double threshold = min_ca + (max_ca - min_ca)*0.8;
                            std::cout << "Ca threshold = " << threshold << "\n";
                            CellProperties calcium_properties(calcium, solution.rGetTimes(),threshold);// Threshold for Ca is tricky...
                            double ca_d50, ca_d90, ca_upstroke, ca_peak;
                            try
                            {
                                ca_d50 = calcium_properties.GetLastActionPotentialDuration(50);
                                ca_d90 = calcium_properties.GetLastActionPotentialDuration(90);
                                ca_upstroke = calcium_properties.GetLastMaxUpstrokeVelocity();
                                ca_peak = calcium_properties.GetLastPeakPotential();
                            }
                            catch (Exception& e)
                            {
                                if (e.GetShortMessage()=="AP did not occur, never exceeded threshold voltage." ||
                                        e.GetShortMessage()=="No full action potential was recorded")
                                {
                                    std::cout << "At Steady pacing of " << hertz << "Hz we did not record any APs\n" << std::flush;
                                    ca_d50 = -1;
                                    ca_d90 = -1;
                                    ca_upstroke = -1;
                                    ca_peak = -1; // Default error/missing data code for matlab analysis.
                                }
                                else
                                {
                                    throw e;
                                }
                            }

                            std::cout << "1Hz Ca Upstroke velocity = " << ca_upstroke << ", Peak Ca = " << ca_peak << ", CaD50 = " << ca_d50 << ", CaD90 = " << ca_d90 << "\n" << std::flush;
                            *steady_calcium_results_file << drug_index << "\t" << drug_conc[conc_index] << "\t" << ca_upstroke << "\t" << ca_peak << "\t" << ca_d50  << "\t" << ca_d90 << "\n";

                            // This uses too much disk space when all models run hERG-only and normal for all drugs.
                            // Create unique filename and write to file... for calcium processing
                            //unsigned output_sig_figs = 6;
                            //std::stringstream filename;
                            //filename << model_name << "_" << drug_index << "_" << drug_conc[conc_index]<< "_"<< j << "_" << hertz;
                            //solution.WriteToFile(steady_state_foldername.str(), filename.str(), "ms", 1, false, output_sig_figs);
                        }
                    } // Stimulus
                } // Steady State

                // Record the state variables at the end of the steady state experiment...
                N_Vector solution_at_1Hz_steady_state = p_model->GetStateVariables(); // This needs changing to a std::vector when using Chaste cells.

                /**
                 * S1-S2 PROTOCOL EXPERIMENT
                 */
                if (s1s2_protocol_experiment)
                {
                    p_model->SetStimulusFunction(p_s1s2_stimulus); // Assign the S1-S2 stimulus to the cell's stimulus

                    // Model is already at the "end of s1 steady state" so commenting out this line.
                    // p_model->Solve(0.0, duration_of_s1, maximum_time_step,printing_time_step );

                    for (unsigned freq=0; freq < p_s1s2_stimulus->GetNumS2FrequencyValues(); freq++)
                    {
                        // Reset the model to the end of S1 (steady state 1Hz pacing).
                        // Take a copy of the saved solution and give it to the cell to play with.
                        p_model->SetStateVariablesUsingACopyOfThisVector(solution_at_1Hz_steady_state);

                        // Update the stimulus to use the next frequency.
                        p_s1s2_stimulus->SetS2ExperimentPeriodIndex(freq);

                        // Solve for two stims with s2 spacing.
                        std::cout << "Solving from " << duration_of_s1 << " to " << duration_of_s1 + 2*s2_periods[freq] << "\n"<<std::flush;
                        OdeSolution solution = p_model->Solve(duration_of_s1, duration_of_s1 + 2*s2_periods[freq], maximum_time_step,printing_time_step );

                        // Get voltage properties
                        std::vector<double> voltages = solution.GetVariableAtIndex(voltage_index); // Voltage should always be index zero
                        CellProperties voltage_properties(voltages, solution.rGetTimes(),voltage_AP_threshold);// Threshold for AP is above -50(ish).
                        std::vector<double> apd50s, apd90s;
                        try
                        {
                            apd50s = voltage_properties.GetAllActionPotentialDurations(50);
                            apd90s = voltage_properties.GetAllActionPotentialDurations(90);
                        }
                        catch (Exception& e)
                        {
                            if (e.GetShortMessage()=="AP did not occur, never exceeded threshold voltage." ||
                                    e.GetShortMessage()=="No full action potential was recorded")
                            {
                                // We can carry on from here OK...
                            }
                            else
                            {
                                throw e;
                            }
                        }

                        if (apd90s.size()==2)
                        {
                            std::cout << "s1s2: Freq = " << s2_periods[freq] << ", APD50 = " << apd50s[0] << ", APD90 = " << apd90s[0] << ", APD50 = " << apd50s[1] << ", APD90 = " << apd90s[1] << "\n" << std::flush;
                            *s1s2_results_file << drug_index << "\t" << drug_conc[conc_index] << "\t" << s2_periods[freq] << "\t" << apd50s[0] << "\t" << apd90s[0] << "\t" << apd50s[1] << "\t" << apd90s[1] << "\n";
                        }
                        else
                        {
                            std::cout << "At S2 PCL = " << s2_periods[freq] << "ms we did not record two independent APs, skipping smaller PCLs\n";
                            break;
                        }

                        // This uses too much disk space when all models run hERG-only and normal for all drugs.
                        // Create unique filename and write to file... for calcium processing
                        //unsigned output_sig_figs = 6;
                        //std::stringstream filename;
                        //filename << model_name << "_" << drug_index << "_" << drug_conc[conc_index]<< "_"<< freq;
                        //solution.WriteToFile(s1s2_foldername.str(), filename.str(), "ms", 1, false, output_sig_figs);
                    } // S2 Freq
                } // S1-S2

                // reset to 1Hz steady state pacing again.
                p_model->SetStateVariablesUsingACopyOfThisVector(solution_at_1Hz_steady_state);

                /**
                * DYNAMIC PACING PROTOCOL EXPERIMENT
                */
                if (dynamic_protocol_experiment)
                {
                    assert(num_pulses>num_pulses_to_analyse);
                    p_model->SetStimulusFunction(p_dynamic_stimulus); // Assign the dynamic stimulus to the cell's stimulus

                    // (We need to keep track of time so that the stimulus works properly.)
                    double run_start_time = 0.0;
                    for (unsigned j=0 ; j<pacing_cycle_lengths.size(); j++)
                    {
                       // Stimulate lots of times to establish the steady state response.
                        p_model->Solve(run_start_time,run_start_time+(num_pulses-num_pulses_to_analyse)*pacing_cycle_lengths[j], maximum_time_step,printing_time_step );
                        run_start_time += (num_pulses-num_pulses_to_analyse)*pacing_cycle_lengths[j]; // Advance time

                        OdeSolution solution = p_model->Solve(run_start_time,run_start_time+num_pulses_to_analyse*pacing_cycle_lengths[j], maximum_time_step,printing_time_step );
                        run_start_time += num_pulses_to_analyse*pacing_cycle_lengths[j]; // Advance time

                        // Get post-processed APDs
                        std::vector<double> voltages = solution.GetVariableAtIndex(voltage_index);
                        CellProperties voltage_properties(voltages, solution.rGetTimes(),voltage_AP_threshold);

                        std::vector<double> apd90s;
                        try
                        {
                            apd90s = voltage_properties.GetAllActionPotentialDurations(90);
                        }
                        catch (Exception& e)
                        {
                            if (e.GetShortMessage()=="AP did not occur, never exceeded threshold voltage." ||
                                    e.GetShortMessage()=="No full action potential was recorded")
                            {
                                std::cout << "At dynamic PCL = " << pacing_cycle_lengths[j] << "ms we did not record any APs, skipping smaller PCLs\n" << std::flush;
                                break;
                            }
                            else
                            {
                                throw e;
                            }
                        }

                        // Check that any action potentials occurred
                        if (apd90s.size()==0)
                        {
                            std::cout << "At dynamic PCL = " << pacing_cycle_lengths[j] << "ms we did not record any APs, skipping smaller PCLs\n" << std::flush;
                            break;
                        }
                        else
                        {
                            std::cout << "PCL = "<< pacing_cycle_lengths[j] << "ms: # action potentials = " << apd90s.size() << "/" << num_pulses_to_analyse << "\n" << std::flush;
                        }

                        // Write out processed results to file.
                        *dynamic_results_file << drug_index << "\t" << drug_conc[conc_index] << "\t" << pacing_cycle_lengths[j];
                        for (unsigned i=0; i< num_pulses_to_analyse; i++)
                        {
                            if (i<apd90s.size())
                            {
                                *dynamic_results_file << "\t" << apd90s[i];
                            }
                            else
                            {
                                *dynamic_results_file << "\t-1";
                            }
                        }
                        *dynamic_results_file << "\n";

                    } // Stimulus
                } // Dynamic

                // Tidy up CVODE vector.
                solution_at_1Hz_steady_state->ops->nvdestroy(solution_at_1Hz_steady_state);

                // Write out results...
                *steady_voltage_results_file << std::flush;
                *steady_calcium_results_file << std::flush;
                *s1s2_results_file << std::flush;
                *dynamic_results_file << std::flush;
            }// Conc
        } // Drug

        // Tidy up
        steady_voltage_results_file->close();
        steady_calcium_results_file->close();
        s1s2_results_file->close();
        dynamic_results_file->close();

        // Clean up memory
        delete p_model;
    }

};


#endif //_TESTMODIFYINGCONDUCTANCES_HPP_

#endif //_CHASTE_CVODE

