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

#ifndef TESTCARDIACSIMULATION_HPP_
#define TESTCARDIACSIMULATION_HPP_

#include <cxxtest/TestSuite.h>

#include "CardiacSimulation.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "CompareHdf5ResultsFiles.hpp"
#include "FileFinder.hpp"

class TestCardiacSimulation : public CxxTest::TestSuite
{
public:

    void TestMono1dSmall() throw(Exception)
    {
        CardiacSimulation simulation("heart/test/data/xml/monodomain1d_small.xml");
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "mono_1d_small", false,
                                                 "SaveMono1D", "SimulationResults", true));
        CardiacSimulation simulation2("heart/test/data/xml/monodomain1d_resume.xml");
    }
    void TestMono2dSmall() throw(Exception)
    {
        CardiacSimulation simulation("heart/test/data/xml/monodomain2d_small.xml");
        CardiacSimulation simulation2("heart/test/data/xml/monodomain2d_resume.xml");
    }
    void TestMono3dSmall() throw(Exception)
    {
        CardiacSimulation simulation("heart/test/data/xml/monodomain3d_small.xml");
        CardiacSimulation simulation2("heart/test/data/xml/monodomain3d_resume.xml");
    }
    void TestBi1dSmall() throw(Exception)
    {
        { CardiacSimulation simulation("heart/test/data/xml/bidomain1d_small.xml"); }
        { CardiacSimulation simulation2("heart/test/data/xml/bidomain1d_resume.xml"); }
        {
            PetscTools::Barrier("TestBi1dSmall-a");
            std::string normal_output_dir = OutputFileHandler::GetChasteTestOutputDirectory();
            setenv("CHASTE_TEST_OUTPUT", (normal_output_dir + "SaveBi1D_checkpoints/0.1ms").c_str(), 1/*Overwrite*/);
            
            try
            {
                // The default resume file specifies a simulation duration of zero.
                // In reality the user should edit the file to specify something sensible...
                TS_ASSERT_THROWS_THIS(CardiacSimulation simulation(OutputFileHandler::GetChasteTestOutputDirectory() + "ResumeParameters.xml"),
                                      "The simulation duration must be positive");
            }
            catch (Exception& e)
            {
                setenv("CHASTE_TEST_OUTPUT", normal_output_dir.c_str(), 1/*Overwrite*/);
                throw e;
            }
            
            PetscTools::Barrier("TestBi1dSmall-b");
            setenv("CHASTE_TEST_OUTPUT", normal_output_dir.c_str(), 1/*Overwrite*/);
        }
    }
    void TestBi2dSmall() throw(Exception)
    {
        CardiacSimulation simulation("heart/test/data/xml/bidomain2d_small.xml");
        CardiacSimulation simulation2("heart/test/data/xml/bidomain2d_resume.xml");
    }
    void TestBi3dSmall() throw(Exception)
    {
        CardiacSimulation simulation("heart/test/data/xml/bidomain3d_small.xml");
        CardiacSimulation simulation2("heart/test/data/xml/bidomain3d_resume.xml");
    }

    void TestCardiacSimulationBasicBidomainShort() throw(Exception)
    {
        // run a bidomain simulation
        CardiacSimulation simulation("heart/test/data/xml/base_bidomain_short.xml");

        // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "base_bidomain_short_results", false,
                   "BaseBidomainShort", "SimulationResults", true));
    }

    void TestCardiacSimulationBasicMonodomainShort() throw(Exception)
    {
        // run a bidomain simulation
        CardiacSimulation simulation("heart/test/data/xml/base_monodomain_short.xml");
        std::string foldername = "BaseMonodomainShort";

       // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "base_monodomain_short_results", false,
                   foldername, "SimulationResults", true));
    }

    void TestCardiacSimulationPostprocessMonodomain() throw(Exception)
    {
        // run a bidomain simulation
        CardiacSimulation simulation("heart/test/data/xml/postprocess_monodomain_short.xml");
        std::string foldername = "PostprocessMonodomainShort";

        // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "postprocess_monodomain_short_results", false,
                   foldername, "SimulationResults", true));
    }

    void TestCardiacSimulationArchiveBidomain() throw(Exception)
    {
        // run a bidomain simulation
        CardiacSimulation simulation("heart/test/data/xml/save_bidomain_short.xml");
        std::string foldername = "SaveBidomainShort";

        // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT(CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "save_bidomain_short_results", false,
                                                foldername, "SimulationResults", true));

        std::string command = "test -e " +  OutputFileHandler::GetChasteTestOutputDirectory() + foldername + "_checkpoints/0.2ms/" + foldername + "_0.2ms/archive.arch.0";
        int return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value,0);
    }

    // requires TestCardiacSimulationArchiveBidomain() to have been run
    void TestCardiacSimulationResumeBidomain() throw(Exception)
    {
        // run a bidomain simulation
        HeartConfig::Instance()->SetSpaceDimension(1);
        CardiacSimulation simulation("heart/test/data/xml/resume_bidomain_short.xml");
        std::string foldername = "SaveBidomainShort";

        // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "resume_bidomain_short_results", false,
                   foldername, "SimulationResults", true));
    }

    void TestCardiacSimulationArchiveMonodomain() throw(Exception)
    {
        // run a bidomain simulation
        CardiacSimulation simulation("heart/test/data/xml/save_monodomain_short.xml");
        std::string foldername = "SaveMonodomainShort";

        // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "save_monodomain_short_results", false,
                   foldername, "SimulationResults", true));

        std::string command = "test -e " +  OutputFileHandler::GetChasteTestOutputDirectory() + foldername + "_checkpoints/0.2ms/" + foldername + "_0.2ms/archive.arch.0";
        int return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value,0);
    }

    // requires TestCardiacSimulationArchiveMonodomain() to have been run
    void TestCardiacSimulationResumeMonodomain() throw(Exception)
    {
        // run a monodomain simulation
        HeartConfig::Instance()->SetSpaceDimension(1);
        CardiacSimulation simulation("heart/test/data/xml/resume_monodomain_short.xml");
        std::string foldername = "SaveMonodomainShort";

        // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "resume_monodomain_short_results", false,
                   foldername, "SimulationResults", true));
    }

    void TestCardiacSimulationArchiveDynamic() throw(Exception)
    {
#ifdef CHASTE_CAN_CHECKPOINT_DLLS
        // run a monodomain simulation
        {
            CardiacSimulation simulation("heart/test/data/xml/save_monodomain_dynamic.xml");
        }
        std::string foldername = "SaveMonodomainDynamic";

        // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "save_monodomain_dynamic", false,
                   foldername, "SimulationResults", true));

        std::string command = "test -e " +  OutputFileHandler::GetChasteTestOutputDirectory() + foldername + "_checkpoints/0.2ms/" + foldername + "_0.2ms/archive.arch.0";
        int return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value,0);

        //resume the simulation
        {
            CardiacSimulation simulation("heart/test/data/xml/resume_monodomain_dynamic.xml");
        }

        // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "resume_monodomain_dynamic", false,
                   foldername, "SimulationResults", true));
#endif // CHASTE_CAN_CHECKPOINT_DLLS
    }

    /**
     * Run TestCardiacSimulationArchiveBidomain on 4 processors to create the archive for this test,
     * and copy it to the repository using:
     *
       scons build=GccOpt_hostconfig,boost=1-33_4 test_suite=heart/test/TestCardiacSimulation.hpp
       cp -r /tmp/$USER/testoutput/SaveBidomainShort_checkpoints/0.2ms heart/test/data/checkpoint_migration_via_xml/
       rm heart/test/data/checkpoint_migration_via_xml/0.2ms/SaveBidomainShort/progress_status.txt
       rm -rf heart/test/data/checkpoint_migration_via_xml/0.2ms/SaveBidomainShort/output
     */
    void TestCardiacSimulationResumeMigration() throw(Exception)
    {
        // We can only load simulations from CHASTE_TEST_OUTPUT, so copy the archives there
        std::string source_directory = "heart/test/data/checkpoint_migration_via_xml/0.2ms/";
        // Clear the target directories
        OutputFileHandler h1("SaveBidomainShort");
        OutputFileHandler h2("SaveBidomainShort_0.2ms");
        if (PetscTools::AmMaster())
        {
            EXPECT0(system, "cp " + source_directory + "SaveBidomainShort/* " + h1.GetOutputDirectoryFullPath());
            EXPECT0(system, "cp " + source_directory + "SaveBidomainShort_0.2ms/* " + h2.GetOutputDirectoryFullPath());
        }
        PetscTools::Barrier("TestCardiacSimulationResumeMigration");

        // Resume
        CardiacSimulation simulation("heart/test/data/xml/resume_migration.xml");
        // Compare results
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "resume_bidomain_short_results", false,
                                                 "SaveBidomainShort", "SimulationResults", true));
    }

    void TestCardiacSimulationPatchwork() throw(Exception)
    {
        OutputFileHandler handler("DynamicallyLoadedModel");
        if (PetscTools::AmMaster())
        {
            // Copy CellML file into output dir
            FileFinder cellml_file("heart/dynamic/luo_rudy_1991.cellml", cp::relative_to_type::chaste_source_root);
            EXPECT0(system, "cp " + cellml_file.GetAbsolutePath() + " " + handler.GetOutputDirectoryFullPath());
        }
        PetscTools::Barrier("TestCardiacSimulationPatchwork");

        CardiacSimulation simulation("heart/test/data/xml/base_monodomain_patchwork.xml");
        std::string foldername = "Patchwork";

        // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT(CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "patchwork_results", false,
                                                foldername, "SimulationResults", true));
    }
    void TestCardiacSimulationKirsten() throw(Exception)
    {
        CardiacSimulation simulation("heart/test/data/xml/base_monodomain_kirsten.xml");
        std::string foldername = "Kirsten";
        TS_ASSERT(CompareFilesViaHdf5DataReaderGlobalNorm("heart/test/data/cardiac_simulations", "Kirsten", false,
                                                foldername, "SimulationResults", true));
    }

    void TestTransmuralCellularheterogeneities() throw(Exception)
    {
        CardiacSimulation simulation("heart/test/data/xml/Transmural_heterogeneities/ChasteParametersCellHeterogeneities.xml");
        std::string foldername = "ChasteResults_heterogeneities";

        TS_ASSERT( CompareFilesViaHdf5DataReaderGlobalNorm("heart/test/data/cardiac_simulations", "transmural_heterogeneities_results", false,
                   foldername, "SimulationResults", true));


    }
    void TestExceptions() throw(Exception)
    {
        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation("heart/test/data/xml/monodomain8d_small.xml"),
                              "Monodomain space dimension not supported: should be 1, 2 or 3");
        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation("heart/test/data/xml/bidomain8d_small.xml"),
                              "Bidomain space dimension not supported: should be 1, 2 or 3");

        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation("heart/test/data/xml/base_monodomain_frankenstein.xml"),
                              "XML parsing error in configuration file: heart/test/data/xml/base_monodomain_frankenstein.xml");

        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation("no file"),
                              "Missing file parsing configuration file: no file");
        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation(""), "No XML file name given");

        FileFinder model("file_does_not_exist.so", cp::relative_to_type::chaste_source_root);
        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation("heart/test/data/xml/missing_dynamic_model.xml"),
                              "Dynamically loadable cell model '" + model.GetAbsolutePath() + "' does not exist.");
#ifndef CHASTE_CAN_CHECKPOINT_DLLS
        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation("heart/test/data/xml/dynamic_checkpoint.xml"),
                              "Checkpointing is not compatible with dynamically loaded cell models on Boost<1.37.");
#endif
    }
};

#endif /*TESTCARDIACSIMULATION_HPP_*/

