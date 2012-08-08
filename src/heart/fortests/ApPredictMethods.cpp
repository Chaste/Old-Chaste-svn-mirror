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

#include "ApPredictMethods.hpp"

#include "AbstractCvodeCell.hpp"
#include "AbstractCardiacCellWithModifiers.hpp"
#include "Exception.hpp"
#include "OutputFileHandler.hpp"
#include "CommandLineArguments.hpp"
#include "FileFinder.hpp"

#include "CellProperties.hpp"
#include "DataStructure.hpp"
#include "AbstractModifier.hpp"
#include "Modifiers.hpp"

#include "ZeroStimulus.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"

#include "Shannon2004Cvode.hpp"
#include "ten_tusscher_model_2006_epiCvode.hpp"
#include "hund_rudy_2004_aCvode.hpp"
#include "mahajan_2008Cvode.hpp"
#include "grandi2010ssCvode.hpp"
#include "ProgressReporter.hpp"

LinearDiscriminantAnalysis ApPredictMethods::LoadLdaFromDrugData()
{
    // We will have to look for the drug data in a couple of places because the executable won't have the Chaste source.
    boost::shared_ptr<FileFinder> p_file_finder = boost::shared_ptr<FileFinder>(new FileFinder("projects/CardiovascRes11/test/drug_data.dat", RelativeTo::ChasteSourceRoot));
    if (!p_file_finder->Exists())
    {
        p_file_finder->SetPath("drug_data.dat", RelativeTo::CWD);
        if (!p_file_finder->Exists())
        {
            EXCEPTION("The file \"drug_data.dat\" should be in the current working directory and is missing.");
        }
    }
    DataStructure drug_data(*p_file_finder);

    matrix<double> cat1and2 (1,1);
    matrix<double> cat3 (1,1);
    matrix<double> cat4 (1,1);
    matrix<double> cat5 (1,1);
    cat1and2(0,0) = cat3(0,0) = cat4(0,0) = cat5(0,0) = DBL_MAX;
    /// \todo THINK ABOUT THIS! Should we add lots of these:
    // cat5(0,0) = 0.0; // Add a control to the classification (no effect = safe)??

    // Generate structures for training data...
    for (unsigned i=0; i<drug_data.GetNumDrugs(); i++)
    {
        double grandi_measure;
        try
        {
            grandi_measure = drug_data.GetGrandiMeasure(i);
        }
        catch (Exception& e)
        {
            assert(e.GetShortMessage()=="No data available on Grandi measure for this drug.");
            continue;
        }

        if (drug_data.GetDrugName(i)=="Ranolazine")
        {   // We aren't currently using Ranolazine for the training.
            continue;
        }

        unsigned redfern = drug_data.GetRedfernCategory(i);
        if (redfern==1 || redfern==2)
        {
            if (cat1and2(0,0) != DBL_MAX)
            {
                cat1and2.resize(cat1and2.size1()+1,1);
            }
            cat1and2(cat1and2.size1()-1,0) = grandi_measure;
        }
        else if (redfern==3)
        {
            if (cat3(0,0) != DBL_MAX)
            {
                cat3.resize(cat3.size1()+1,1);
            }
            cat3(cat3.size1()-1,0) = grandi_measure;
        }
        else if (redfern==4)
        {
            if (cat4(0,0) != DBL_MAX)
            {
                cat4.resize(cat4.size1()+1,1);
            }
            cat4(cat4.size1()-1,0) = grandi_measure;
        }
        else
        {
            assert(redfern==5);
            if (cat5(0,0) != DBL_MAX)
            {
                cat5.resize(cat5.size1()+1,1);
            }
            cat5(cat5.size1()-1,0) = grandi_measure;
        }
    }

    std::vector<matrix<double> > training;
    training.push_back(cat1and2);
    training.push_back(cat3);
    training.push_back(cat4);
    training.push_back(cat5);

    LinearDiscriminantAnalysis lda(training);
    return lda;
}


void ApPredictMethods::MakeTorsadePredictions()
{
    if (mAPD90s.size()==0)
    {
        EXCEPTION("APDs do not appear to have been recorded.");
    }
    mTorsadePredictions.reserve(mAPD90s.size());

    // Work out the Measure
    std::vector<double> largest_percent_change;
    largest_percent_change.push_back(0.0); // First is always control
    for (unsigned i=1; i<mAPD90s.size(); i++)
    {
        double percent = 100*(mAPD90s[i]-mAPD90s[0])/mAPD90s[0];
        largest_percent_change.push_back(percent);
    }

    // Second pass to check for largest +ve effect at low dose...
    // This makes the measure "Largest Effect Dose" rather than just "dose".
    for (unsigned i=1; i<mAPD90s.size(); i++)
    {
        // Check over lower doses
        for (unsigned j=0; j<i; j++)
        {
            // If the lower dose made the AP longer than control (and longer than this dose) then replace this APD with that one.
            if (largest_percent_change[j] > 0 && largest_percent_change[j] > largest_percent_change[i])
            {
                largest_percent_change[i] = largest_percent_change[j];
            }
        }
    }

    // Perform the LDA
    LinearDiscriminantAnalysis lda = LoadLdaFromDrugData();
    for (unsigned i=0; i<mAPD90s.size(); i++)
    {
        vector<double> test_point(1);
        test_point(0) = largest_percent_change[i];
        mTorsadePredictions.push_back(lda.ClassifyThisPoint(test_point)+2u); // We add two because our redfern categories start at 2 not 0.
    }
}


ApPredictMethods::ApPredictMethods(bool torsadePredict)
{
    mAPD90s.clear();
    CommandLineArguments* p_args = CommandLineArguments::Instance();
    char **argv = *(p_args->p_argv); // is a char** of them.
    unsigned model_index;
    double hERG_IC50, Na_IC50, CaL_IC50, top_dose;
    std::string program_name;
    std::string output_folder_base;
    if (torsadePredict)
    {
        model_index = 5u;// Hardcoded to Grandi model.
        hERG_IC50 = atof(argv[1]);
        Na_IC50 = atof(argv[2]);
        CaL_IC50 = atof(argv[3]);
        top_dose = atof(argv[4]);
        program_name = "Torsade PreDiCT";
        output_folder_base = "TorsadePredict_output";
    }
    else
    {
        model_index = atoi(argv[1]);
        hERG_IC50 = atof(argv[2]);
        Na_IC50 = atof(argv[3]);
        CaL_IC50 = atof(argv[4]);
        top_dose = atof(argv[5]);
        program_name = "Action Potential PreDiCT";
        output_folder_base = "ApPredict_output/";
        std::cout << "* model_index = " << model_index << "\n";
    }
    std::cout << "* hERG IC50   = " << hERG_IC50 << " nM\n"
                 "* Na IC50     = " << Na_IC50 << " nM\n"
                 "* CaL IC50    = " << CaL_IC50 << " nM\n"
                 "* max free plasma concentration = " << top_dose << " nM\n" << std::flush;

    double s_start = 5.0;           // ms
    double printing_time_step = 1.0;  // ms

    // PARAMETERS FOR THE STEADY STATE EXPERIMENTS
    // This is a regular stimulus for steady-state pacing experiments
    boost::shared_ptr<ZeroStimulus> p_zero_stimulus(new ZeroStimulus());
    const unsigned num_stims = 1000u;

    /// Cvode cells use a CVODE solver regardless of which standard solver is passed in.
    boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

    // Pick the model...
    // 1= Shannon, 2=TenTusscher, 3 = Mahajan, 4 = HundRudy, 5=Grandi
    AbstractCardiacCellWithModifiers<AbstractCvodeCell>* p_model;

    switch (model_index)
    {
        case 1u:
            p_model = new CellShannon2004FromCellMLCvode(p_solver, p_zero_stimulus);
            break;
        case 2u:
            p_model = new Cellten_tusscher_model_2006_epiFromCellMLCvode(p_solver, p_zero_stimulus);
            break;
        case 3u:
            p_model = new Cellmahajan_2008FromCellMLCvode(p_solver, p_zero_stimulus);
            break;
        case 4u:
            p_model = new Cellhund_rudy_2004_aFromCellMLCvode(p_solver, p_zero_stimulus);
            break;
        case 5u:
            p_model = new Cellgrandi2010ssFromCellMLCvode(p_solver, p_zero_stimulus);
            break;
        default:
            EXCEPTION("No model matches this index");
    }

    // Tell the model to use the default stimulus amplitude/duration/period assigned to it in CellML.
    p_model->UseCellMLDefaultStimulus();
    boost::shared_ptr<RegularStimulus> p_reg_stim =
                 boost::static_pointer_cast<RegularStimulus>(p_model->GetStimulusFunction());
    double s1_period = p_reg_stim->GetPeriod();
    double s_duration = p_reg_stim->GetDuration();
    double hertz = 1000.0/s1_period;
    p_reg_stim->SetStartTime(s_start);
    double maximum_time_step; // ms
    if (printing_time_step > s_duration)
    {
        maximum_time_step = s_duration;
    }
    else
    {
        maximum_time_step = printing_time_step;
    }

    boost::shared_ptr<const AbstractOdeSystemInformation> p_ode_info = p_model->GetSystemInformation();
    std::string model_name = p_model->GetSystemName();

    // Set up foldernames for each model and protocol set.
    std::string foldername = output_folder_base;

    // Make and clean the above directories.
    OutputFileHandler file_handler(foldername);

    unsigned num_doses = 11;
    ProgressReporter progress_reporter(foldername, 0.0, (double)(num_doses));
    progress_reporter.PrintInitialising();

    // Open files and write headers
    out_stream steady_voltage_results_file_html =  file_handler.OpenOutputFile("voltage_results.html");

    out_stream steady_voltage_results_file =  file_handler.OpenOutputFile("voltage_results.dat");
    *steady_voltage_results_file << "Concentration(nM)\tUpstrokeVelocity(mV/ms)\tPeakVm(mV)\tAPD50(ms)\tAPD90(ms)\n";
    *steady_voltage_results_file_html << "<html>\n<head><title>" << program_name << " results</title></head>\n";
    *steady_voltage_results_file_html << "<STYLE TYPE=\"text/css\">\n<!--\nTD{font-family: Arial; font-size: 12pt;}\n--->\n</STYLE>\n";
    *steady_voltage_results_file_html << "<body>\n";
    *steady_voltage_results_file_html << "<table width=\"60%\" style=\"background-color:white\" border=\"1\" cellpadding=\"2\" cellspacing=\"0\">\n";
    *steady_voltage_results_file_html << "<tr><td>Concentration (nM)</td><td>Upstroke Velocity (mV/ms)</td><td>Peak membrane voltage (mV)</td><td>APD50 (ms)</td><td>APD90 (ms)</td></tr>\n"; // Header line


    unsigned voltage_index = p_ode_info->GetStateVariableIndex("membrane_voltage");

    /**
     * DECIDE WHAT CONCENTRATIONS TO TEST AT
     */
    for (unsigned i=0 ; i<num_doses ;i++)
    {
        mConcs.push_back(((double)(i)/(double)(num_doses-1))*top_dose);
    }

    /**
     * START LOOP OVER EACH CONCENTRATION TO TEST WITH
     */
    for (unsigned conc_index=0; conc_index<mConcs.size(); conc_index++)
    {
        progress_reporter.Update((double)(conc_index));
        std::cout << "Drug Conc = " << mConcs[conc_index] << " nM\n" << std::flush;

        // Here we calculate the proportion of the different channels which are still active
        // (at this concentration of this drug)
        double gNa_factor = DataStructure::CalculateConductanceFactor(mConcs[conc_index],Na_IC50);
        double gCaL_factor = DataStructure::CalculateConductanceFactor(mConcs[conc_index],CaL_IC50);
        double gKr_factor = DataStructure::CalculateConductanceFactor(mConcs[conc_index],hERG_IC50);

        std::cout << "gNa factor = " << gNa_factor << "\n" << std::flush;
        std::cout << "gCaL factor = " << gCaL_factor << "\n" << std::flush;
        std::cout << "gKr factor = " << gKr_factor << "\n" << std::flush;

        // 'Factor Modifier' classes allow us to intervene in the right hand side of the ODE,
        // multiplying the usual values of the labelled parameter by a scaling factor.
        boost::shared_ptr<AbstractModifier> a = boost::shared_ptr<AbstractModifier>(new FactorModifier(gNa_factor));
        boost::shared_ptr<AbstractModifier> b = boost::shared_ptr<AbstractModifier>(new FactorModifier(gCaL_factor));
        boost::shared_ptr<AbstractModifier> c = boost::shared_ptr<AbstractModifier>(new FactorModifier(gKr_factor));

        // The following names are fixed and correspond to metadata tags.
        p_model->SetModifier("membrane_fast_sodium_current_conductance", a);
        p_model->SetModifier("membrane_L_type_calcium_current_conductance", b);
        p_model->SetModifier("membrane_rapid_delayed_rectifier_potassium_current_conductance", c);

        /**
         * STEADY STATE PACING EXPERIMENT
         *
         * Stimulate many times to establish the steady state response.
         * Do most of the calculations in one go to avoid overhead setting up Cvode...
         */
        p_model->SetMaxSteps(1e7);
        p_model->Solve(0.0, s1_period*(double)(num_stims), maximum_time_step);

        // Set the max time step to be the printing time step (1ms).
        // Do one final pace that we will analyse...
        OdeSolution solution = p_model->Solve(s1_period*(double)num_stims, s1_period*(double)(num_stims+1), maximum_time_step, printing_time_step);

        // Get voltage properties
        std::vector<double> voltages = solution.GetVariableAtIndex(voltage_index); // Voltage should always be zero
        CellProperties voltage_properties(voltages, solution.rGetTimes(),-50);// Threshold for AP is above -50ish.

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
        *steady_voltage_results_file_html << "<tr><td>"<< mConcs[conc_index] << "</td><td>" << upstroke << "</td><td>" << peak << "</td><td>" << apd50 << "</td><td>" << apd90 << "</td></tr>\n";
        *steady_voltage_results_file << mConcs[conc_index] << "\t" << upstroke << "\t" << peak << "\t" << apd50  << "\t" << apd90 << "\n";

        if (torsadePredict) // If we are going to make Torsade predictions then assign APD90 to a member variable
        {
            mAPD90s.push_back(apd90);
        }

        // Create unique filename and write to file...
        std::stringstream filename;
        filename << "conc_" << conc_index << "_voltage_trace.dat";
        out_stream voltage_results_file =  file_handler.OpenOutputFile(filename.str());
        *voltage_results_file << "Time(ms)\tMembrane_Voltage(mV)\n";
        double last_voltage_printed;
        bool printed_last = true; // This is to make sure large jumps print out the step before.
        for (unsigned i=0; i<voltages.size(); i++)
        {
            // A new bit of code to downsample the output so flot doesn't timeout loading large datasets
            if (i>0 && i<voltages.size()-1) // if we aren't at the beginning or the end of the trace
            {
                if (fabs(voltages[i] - last_voltage_printed) < 1.0)
                {   // Only output a point if the voltage has changed by 1 mV.
                    printed_last = false;
                    continue;
                }
            }
            if(!printed_last) // We want the point before printing too to avoid large interpolations.
            {
                *voltage_results_file << solution.rGetTimes()[i-1] - solution.rGetTimes()[0] - s_start << "\t" << voltages[i-1] << "\n";
            }
            *voltage_results_file << solution.rGetTimes()[i] - solution.rGetTimes()[0] - s_start << "\t" << voltages[i] << "\n";
            printed_last = true;
            last_voltage_printed = voltages[i];
        }
        voltage_results_file->close();
        *steady_voltage_results_file << std::flush;
    }// Conc

    // Tidy up
    progress_reporter.PrintFinalising();
    *steady_voltage_results_file_html << "</table>\n</body>\n</html>\n";
    steady_voltage_results_file_html->close();
    steady_voltage_results_file->close();

    if (torsadePredict)
    {
        MakeTorsadePredictions();

        // Open an output file for the Torsade results
        OutputFileHandler tdp_handler(output_folder_base, false); // don't wipe the directory!
        out_stream torsade_results_file = file_handler.OpenOutputFile("tdp_results.html");
        *torsade_results_file << "<html>\n<head><title>Torsade preDiCT Results</title></head>\n";
        *torsade_results_file << "<STYLE TYPE=\"text/css\">\n<!--\nTD{font-family: Arial; font-size: 12pt;}\n--->\n</STYLE>\n";
        *torsade_results_file << "<body>\n";
        *torsade_results_file << "<table width=\"60%\" style=\"background-color:white\" border=\"1\" cellpadding=\"2\" cellspacing=\"0\">\n";
        *torsade_results_file << "<tr><td>Concentration (nM)</td><td>APD90 (ms)</td><td>Risk Category Prediction</td></tr>\n"; // Header line

        for (unsigned i=0; i<mTorsadePredictions.size(); i++)
        {
            std::cout << "Conc = " << mConcs[i] << "nM, APD90 = " << mAPD90s[i] << "ms, risk prediction =  " << mTorsadePredictions[i] << "\n" << std::flush;
            *torsade_results_file << "<tr><td>"<< mConcs[i] << "</td><td>" << mAPD90s[i] << "</td><td>" << mTorsadePredictions[i] << "</td></tr>\n";
        }
        *torsade_results_file << "</table>\n</body>\n</html>\n";
        torsade_results_file->close();
    }
    // Clean up memory
    delete p_model;
}
