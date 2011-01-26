/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef TESTPYCMLNIGHTLY_HPP_
#define TESTPYCMLNIGHTLY_HPP_

#include <cxxtest/TestSuite.h>

#include <iostream>

#include "AbstractCardiacCell.hpp"
#include "AbstractCardiacCellInterface.hpp"
#include "AbstractCvodeCell.hpp"

#include "DynamicLoadingHelperFunctions.hpp"

#include "DynamicCellModelLoader.hpp"
#include "DynamicModelLoaderRegistry.hpp"
#include "FileFinder.hpp"
#include "CellMLToSharedLibraryConverter.hpp"

#include "CellProperties.hpp"
#include "HeartConfig.hpp"
#include "Warnings.hpp"
#include "RunAndCheckIonicModels.hpp"

#include "PetscSetupAndFinalize.hpp"

/**
 * Test PyCml functionality by dynamically loading (and hence converting) a wide
 * range of cell models.
 *
 * May need a test-suite setup or similar to define model-specific parameters?
 * Should we pick up the list of models by reading the folder heart/test/data/cellml?
 */
class TestPyCmlNightly : public CxxTest::TestSuite
{
private:
    template<typename VECTOR>
    double GetAttribute(AbstractParameterisedSystem<VECTOR>* pSystem,
                        const std::string& rAttrName,
                        double defaultValue)
    {
        assert(pSystem);
        double attr_value;
        if (pSystem->HasAttribute(rAttrName))
        {
            attr_value = pSystem->GetAttribute(rAttrName);
        }
        else
        {
            attr_value = defaultValue;
        }
        return attr_value;
    }

    double GetAttribute(boost::shared_ptr<AbstractCardiacCellInterface> pCell,
                        const std::string& rAttrName,
                        double defaultValue)
    {
        AbstractParameterisedSystem<std::vector<double> >* p_normal
            = dynamic_cast<AbstractParameterisedSystem<std::vector<double> >*>(pCell.get());
        if (p_normal)
        {
            return GetAttribute(p_normal, rAttrName, defaultValue);
        }
#ifdef CHASTE_CVODE
        else
        {
            AbstractParameterisedSystem<N_Vector>* p_cvode
                = dynamic_cast<AbstractParameterisedSystem<N_Vector>*>(pCell.get());
            return GetAttribute(p_cvode, rAttrName, defaultValue);
        }
#endif
        EXCEPTION("Unrecognised cell type.");
    }

    void Simulate(const std::string& rOutputDirName,
                  const std::string& rModelName,
                  boost::shared_ptr<AbstractCardiacCellInterface> pCell)
    {
        double end_time = GetAttribute(pCell, "SuggestedCycleLength", 700.0); // ms
        if (pCell->GetSolver())
        {
            double dt = GetAttribute(pCell, "SuggestedForwardEulerTimestep", 0.0);
            if (dt > 0.0)
            {
                pCell->SetTimestep(dt);
            }
        }
        double sampling_interval = 2.0; // ms
        OdeSolution solution = pCell->Compute(0.0, end_time, sampling_interval);
        const unsigned output_freq = 5; // Only output every N samples
        solution.WriteToFile(rOutputDirName, rModelName, "ms", output_freq, false);
        // Check an AP was produced
        std::vector<double> voltages = solution.GetVariableAtIndex(pCell->GetVoltageIndex());
        CellProperties props(voltages, solution.rGetTimes());
        props.GetLastActionPotentialDuration(90.0); // Don't catch the exception here if it's thrown
        // Compare against saved results
        CheckResults(rModelName, voltages, solution.rGetTimes(), output_freq);
    }

    void CheckResults(const std::string& rModelName,
                      std::vector<double>& rVoltages,
                      std::vector<double>& rTimes,
                      unsigned outputFreq,
                      double tolerance=1.0)
    {
        // Read data entries for the reference file
        ColumnDataReader data_reader("heart/test/data/cellml", rModelName, false);
        std::vector<double> valid_times = data_reader.GetValues("Time");
        std::vector<double> valid_voltages = GetVoltages(data_reader);

        TS_ASSERT_EQUALS(rTimes.size(), (valid_times.size()-1)*outputFreq+1);
        for (unsigned i=0; i<valid_times.size(); i++)
        {
            TS_ASSERT_DELTA(rTimes[i*outputFreq], valid_times[i], 1e-12);
            TS_ASSERT_DELTA(rVoltages[i*outputFreq], valid_voltages[i], tolerance);
        }
    }

    void RunTests(const std::string& rOutputDirName,
                  const std::vector<std::string>& rModels,
                  const std::vector<std::string>& rArgs,
                  bool testLookupTables=false,
                  double tableTestV=-1000)
    {
        OutputFileHandler handler(rOutputDirName); // Clear folder
        std::vector<std::string> failures;
        for (unsigned i=0; i<rModels.size(); ++i)
        {
            try
            {
                RunTest(rOutputDirName, rModels[i], rArgs, testLookupTables, tableTestV);
            }
            catch (const Exception& e)
            {
                failures.push_back(rModels[i]);
                TS_FAIL("Failure testing cell model " + rModels[i] + ": " + e.GetMessage());
            }
        }

        if (!failures.empty())
        {
            std::cout << failures.size() << " models failed for " << rOutputDirName << ":" << std::endl;
            for (unsigned i=0; i<failures.size(); ++i)
            {
                std::cout << "   " << failures[i] << std::endl;
            }
        }
    }

    void RunTest(const std::string& rOutputDirName,
                 const std::string& rModelName,
                 const std::vector<std::string>& rArgs,
                 bool testLookupTables=false,
                 double tableTestV=-1000)
    {
        // Copy CellML file (and .out if present) into output dir
        OutputFileHandler handler(rOutputDirName, false);
        FileFinder cellml_file("heart/test/data/cellml/" + rModelName + ".cellml", RelativeTo::ChasteSourceRoot);
        CopyFile(handler, cellml_file);
        FileFinder out_file("heart/test/data/cellml/" + rModelName + ".out", RelativeTo::ChasteSourceRoot);
        if (out_file.Exists())
        {
            CopyFile(handler, out_file);
        }

        // Create options file
        if (!rArgs.empty())
        {
            CreateOptionsFile(handler, rModelName, rArgs);
        }

        // Do the conversion
        CellMLToSharedLibraryConverter converter;
        FileFinder copied_file(rOutputDirName + "/" + rModelName + ".cellml", RelativeTo::ChasteTestOutput);
        DynamicCellModelLoader* p_loader = converter.Convert(copied_file);
        // Apply a stimulus of -40 uA/cm^2 - should work for all models
        boost::shared_ptr<AbstractCardiacCellInterface> p_cell(CreateCellWithStandardStimulus(*p_loader, -40.0));

        // Check lookup tables exist if they should
        if (testLookupTables && rModelName != "hodgkin_huxley_squid_axon_model_1952_modified")
        {
            double v = p_cell->GetVoltage();
            p_cell->SetVoltage(tableTestV);
            TS_ASSERT_THROWS_CONTAINS(p_cell->GetIIonic(), "membrane_voltage outside lookup table range");
            p_cell->SetVoltage(v);
        }
        Simulate(rOutputDirName, rModelName, p_cell);
        Warnings::NoisyDestroy(); // Print out any warnings now, not at program exit
    }

    void AddAllModels(std::vector<std::string>& rModels)
    {
        rModels.push_back("aslanidi_model_2009");
        rModels.push_back("beeler_reuter_model_1977");
        rModels.push_back("bondarenko_model_2004_apex");
        rModels.push_back("courtemanche_ramirez_nattel_model_1998");
        rModels.push_back("decker_2009");
        rModels.push_back("demir_model_1994");
        rModels.push_back("dokos_model_1996");
        rModels.push_back("earm_noble_model_1990");
        rModels.push_back("espinosa_model_1998_normal");
        rModels.push_back("fink_noble_giles_model_2008");
        rModels.push_back("grandi2010ss");
        rModels.push_back("hilgemann_noble_model_1987");
        rModels.push_back("hodgkin_huxley_squid_axon_model_1952_modified");
        rModels.push_back("hund_rudy_2004_a");
        rModels.push_back("iribe_model_2006_without_otherwise_section");
        rModels.push_back("iyer_model_2004");
        rModels.push_back("iyer_model_2007");
        rModels.push_back("jafri_rice_winslow_model_1998");
        rModels.push_back("kurata_model_2002");
        rModels.push_back("livshitz_rudy_2007");
        rModels.push_back("mahajan_2008");
        rModels.push_back("matsuoka_model_2003");
        rModels.push_back("noble_model_1991");
        rModels.push_back("noble_model_1998");
        rModels.push_back("noble_noble_SAN_model_1984");
        rModels.push_back("noble_SAN_model_1989");
        rModels.push_back("nygren_atrial_model_1998");
        rModels.push_back("pandit_model_2001_epi");
        rModels.push_back("priebe_beuckelmann_model_1998");
        rModels.push_back("sakmann_model_2000_epi");
        rModels.push_back("stewart_zhang_model_2008_ss");
        rModels.push_back("ten_tusscher_model_2004_endo");
        rModels.push_back("ten_tusscher_model_2004_epi");
        rModels.push_back("ten_tusscher_model_2006_epi");
        rModels.push_back("viswanathan_model_1999_epi");
        rModels.push_back("winslow_model_1999");
        rModels.push_back("zhang_SAN_model_2000_0D_capable");
        rModels.push_back("zhang_SAN_model_2000_all");
    }

public:
    void TestNormalCells() throw (Exception)
    {
        std::cout << "Search for ': ***', 'Error', or 'Failed' to find problems." << std::endl;

        std::string dirname("TestPyCmlNightlyNormal");
        std::vector<std::string> args;
        args.push_back("--Wu");
        std::vector<std::string> models;
        AddAllModels(models);
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.005, 0.1, 1.0);
        RunTests(dirname, models, args);
    }

    void TestOptimisedCells() throw (Exception)
    {
        std::string dirname("TestPyCmlNightlyOpt");
        std::vector<std::string> args;
        args.push_back("--Wu");
        args.push_back("--opt");
        std::vector<std::string> models;
        AddAllModels(models);
        RunTests(dirname, models, args, true);
    }

    void TestCvodeCells() throw (Exception)
    {
        std::string dirname("TestPyCmlNightlyCvode");
        std::vector<std::string> args;
        args.push_back("--Wu");
        args.push_back("--cvode");
        std::vector<std::string> models;
        AddAllModels(models);
        RunTests(dirname, models, args);
    }

    void TestBackwardEulerCells() throw (Exception)
    {
        std::string dirname("TestPyCmlNightlyBE");
        std::vector<std::string> args;
        args.push_back("--Wu");
        args.push_back("--backward-euler");
        std::vector<std::string> models;

        models.push_back("aslanidi_model_2009");
        models.push_back("beeler_reuter_model_1977");
        models.push_back("courtemanche_ramirez_nattel_model_1998");
        models.push_back("demir_model_1994");
        models.push_back("dokos_model_1996");
        models.push_back("earm_noble_model_1990");
        models.push_back("espinosa_model_1998_normal");
        models.push_back("fink_noble_giles_model_2008");
//        models.push_back("grandi2010ss");
        models.push_back("hodgkin_huxley_squid_axon_model_1952_modified");
        models.push_back("kurata_model_2002");
        models.push_back("livshitz_rudy_2007");
        models.push_back("matsuoka_model_2003");
        models.push_back("noble_model_1991");
        models.push_back("noble_model_1998");
        models.push_back("noble_noble_SAN_model_1984");
        models.push_back("noble_SAN_model_1989");
        models.push_back("nygren_atrial_model_1998");
        models.push_back("sakmann_model_2000_epi");
        models.push_back("ten_tusscher_model_2006_epi");
        models.push_back("zhang_SAN_model_2000_0D_capable");
        models.push_back("zhang_SAN_model_2000_all");

        /* The following models contained 'diff' in the .out file:
            decker_2009
            hilgemann_noble_model_1987
            hund_rudy_2004_a
            mahajan_2008
            priebe_beuckelmann_model_1998
            ten_tusscher_model_2004_endo
            ten_tusscher_model_2004_epi
            viswanathan_model_1999_epi
            winslow_model_1999
         */

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.1, 1.0);
        RunTests(dirname, models, args, true);

        dirname = dirname + "-difficult";
        models.clear();
        models.push_back("bondarenko_model_2004_apex");
        models.push_back("iyer_model_2004");
        models.push_back("iyer_model_2007");
        models.push_back("jafri_rice_winslow_model_1998");
        models.push_back("pandit_model_2001_epi");
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.001, 0.1, 1.0);
        RunTests(dirname, models, args, true);
    }
};

#endif // TESTPYCMLNIGHTLY_HPP_
