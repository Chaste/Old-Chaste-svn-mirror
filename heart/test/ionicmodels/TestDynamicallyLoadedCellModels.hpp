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

#include "ChasteSerialization.hpp"
#ifdef CHASTE_CAN_CHECKPOINT_DLLS
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // CHASTE_CAN_CHECKPOINT_DLLS

#include "RunAndCheckIonicModels.hpp"
#include "SimpleStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "DynamicCellModelLoader.hpp"
#include "DynamicModelLoaderRegistry.hpp"
#include "ChasteBuildRoot.hpp"
#include "HeartConfig.hpp"
#include "FileFinder.hpp"
#include "HeartFileFinder.hpp"
#include "CellMLToSharedLibraryConverter.hpp"
#include "AbstractDynamicallyLoadableEntity.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestDynamicallyLoadedCellModels : public CxxTest::TestSuite
{
private:

    void RunLr91Test(DynamicCellModelLoader& rLoader,
                     unsigned vIndex=4u)
    {
        AbstractCardiacCell* p_cell = CreateLr91CellFromLoader(rLoader, vIndex);
        SimulateLr91AndCompare(p_cell);
        delete p_cell;
    }

    void SimulateLr91AndCompare(AbstractCardiacCell* pCell)
    {
        double end_time = 1000.0; //One second in milliseconds
        // Solve and write to file
        clock_t ck_start = clock();
        RunOdeSolverWithIonicModel(pCell,
                                   end_time,
                                   "DynamicallyLoadableLr91");
        clock_t ck_end = clock();
        double forward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        std::cout << "\n\tSolve time: " << forward << std::endl;

        // Compare with 'normal' LR91 model results
        CheckCellModelResults("DynamicallyLoadableLr91", "Lr91DelayedStim");

        // Test GetIIonic against hardcoded result from TestIonicModels.hpp
        RunOdeSolverWithIonicModel(pCell,
                                   60.0,
                                   "DynamicallyLoadableLr91GetIIonic");
        TS_ASSERT_DELTA(pCell->GetIIonic(), 1.9411, 1e-3);
    }
    
    AbstractCardiacCell* CreateLr91CellFromLoader(DynamicCellModelLoader& rLoader,
                                                  unsigned vIndex=4u)
    {
        // Set stimulus
        double magnitude = -25.5;
        double duration  = 2.0 ;  // ms
        double when = 50.0; // ms

        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude, duration, when));
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        // Load the cell model dynamically
        AbstractCardiacCell* p_cell = rLoader.CreateCell(p_solver, p_stimulus);

        // Simple sanity checks
        TS_ASSERT_EQUALS(p_cell->GetVoltageIndex(), vIndex);
        AbstractDynamicallyLoadableEntity* p_entity = dynamic_cast<AbstractDynamicallyLoadableEntity*>(p_cell);
        if (p_entity != NULL)
        {
            TS_ASSERT_EQUALS(&rLoader, p_entity->GetLoader());
            std::cout << "Cell inherits from AbstractDynamicallyLoadableEntity" << std::endl;
        }
        
        return p_cell;
    }
    
public:
    /**
     * This is based on TestOdeSolverForLR91WithDelayedSimpleStimulus from
     * TestIonicModels.hpp.
     */
    void TestDynamicallyLoadedLr91() throw(Exception)
    {
        // Load the cell model dynamically
        std::string model_name = "libDynamicallyLoadableLr91.so";
        // All tests use the registry, as not doing so can lead to segfaults...
        DynamicCellModelLoader* p_loader = DynamicModelLoaderRegistry::Instance()->GetLoader(
            ChasteComponentBuildDir("heart") + "dynamic/" + model_name);
        RunLr91Test(*p_loader);

        // The .so also gets copied into the source folder
        DynamicCellModelLoader* p_loader2 = DynamicModelLoaderRegistry::Instance()->GetLoader(
            std::string(ChasteBuildRootDir()) + "heart/dynamic/" + model_name);
        RunLr91Test(*p_loader2);
    }

    /**
     * Currently the build system will not automatically regenerate the C++ code from CellML, so you
     * need to do that step manually:
     *     cdchaste
     *     ./python/ConvertCellModel.py -A -y heart/dynamic/luo_rudy_1991.cellml
     *
     * \todo #1164 Remove this test?  TestCellmlConverter makes it slightly redundant,
     * apart from the optimised version.
     */
    void TestLr91FromCellML() throw(Exception)
    {
        FileFinder model("heart/dynamic/libluo_rudy_1991.so", RelativeTo::ChasteSourceRoot);
        DynamicCellModelLoader* p_loader = DynamicModelLoaderRegistry::Instance()->GetLoader(model);
        RunLr91Test(*p_loader, 0u);

        FileFinder model_opt("heart/dynamic/libluo_rudy_1991Opt.so", RelativeTo::ChasteSourceRoot);
        DynamicCellModelLoader* p_loader_opt = DynamicModelLoaderRegistry::Instance()->GetLoader(model_opt);
        RunLr91Test(*p_loader_opt, 0u);
        
        // Coverage
        TS_ASSERT_EQUALS(p_loader->GetLoadableModulePath(), model.GetAbsolutePath());
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

    void TestLoadingViaXml() throw(Exception)
    {
        // Fake content from an XML file.
        std::string model_name = "heart/dynamic/libDynamicallyLoadableLr91.so";
        cp::path_type so_path(model_name);
        so_path.relative_to(cp::relative_to_type::chaste_source_root);
        cp::dynamically_loaded_ionic_model_type dynamic_elt(so_path);
        cp::ionic_model_selection_type ionic_model;
        ionic_model.Dynamic(dynamic_elt);

        // Now mock up what HeartConfigRelatedCellFactory does
        TS_ASSERT(ionic_model.Dynamic().present());
        if (ionic_model.Dynamic().present())
        {
            HeartFileFinder file_finder(ionic_model.Dynamic()->Path());
            TS_ASSERT(file_finder.Exists());
            DynamicCellModelLoader* p_loader = DynamicModelLoaderRegistry::Instance()->GetLoader(file_finder);
            RunLr91Test(*p_loader);
        }
    }

    void TestCellmlConverter() throw(Exception)
    {
        // Copy CellML file into output dir
        std::string dirname = "TestCellmlConverter";
        OutputFileHandler handler(dirname);
        if (PetscTools::AmMaster())
        {
            FileFinder cellml_file("heart/dynamic/luo_rudy_1991.cellml", RelativeTo::ChasteSourceRoot);
            EXPECT0(system, "cp " + cellml_file.GetAbsolutePath() + " " + handler.GetOutputDirectoryFullPath());
        }
        PetscTools::Barrier("TestCellmlConverter_cp");

        CellMLToSharedLibraryConverter converter;

        // Convert a real CellML file
        FileFinder cellml_file(dirname + "/luo_rudy_1991.cellml", RelativeTo::ChasteTestOutput);
        TS_ASSERT(cellml_file.Exists());
        FileFinder so_file(dirname + "/libluo_rudy_1991.so", RelativeTo::ChasteTestOutput);
        TS_ASSERT(!so_file.Exists());
        DynamicCellModelLoader* p_loader = converter.Convert(cellml_file);
        TS_ASSERT(so_file.Exists());
        TS_ASSERT(so_file.IsNewerThan(cellml_file));
        // Converting a .so should be a "no-op"
        DynamicCellModelLoader* p_loader2 = converter.Convert(so_file);
        TS_ASSERT(so_file.Exists());
        TS_ASSERT(p_loader2 == p_loader);
        RunLr91Test(*p_loader, 0u);

        // Cover exceptions
        std::string file_name = "test";
        FileFinder no_ext(dirname + "/" + file_name, RelativeTo::ChasteTestOutput);
        TS_ASSERT_THROWS_THIS(converter.Convert(no_ext), "Dynamically loadable cell model '"
                              + no_ext.GetAbsolutePath() + "' does not exist.");

        PetscTools::Barrier("TestCellmlConverter_pre_touch");
        if (PetscTools::AmMaster())
        {
            out_stream fp = handler.OpenOutputFile(file_name);
            fp->close();
        }
        PetscTools::Barrier("TestCellmlConverter_post_touch");

        TS_ASSERT_THROWS_THIS(converter.Convert(no_ext), "File does not have an extension: " + no_ext.GetAbsolutePath());
        FileFinder unsupp_ext("global/src/FileFinder.hpp", RelativeTo::ChasteSourceRoot);
        TS_ASSERT_THROWS_THIS(converter.Convert(unsupp_ext), "Unsupported extension '.hpp' of file '"
                              + unsupp_ext.GetAbsolutePath() + "'; must be .so or .cellml");

        EXPECT0(chdir, "heart"); // The ConvertCellModel.py script in ConvertCellmlToSo() should only work from chaste source directory.
        if (PetscTools::AmMaster())
        {
            // Having the 'rm' after the 'chdir' cunningly double-checks that the FileFinder has really given us an absolute path
            EXPECT0(system, "rm " + so_file.GetAbsolutePath()); // Make sure the conversion is re-run
        }
        PetscTools::Barrier("TestCellmlConverter_rm");
        TS_ASSERT_THROWS_THIS(converter.Convert(cellml_file, false),
                              "Unable to convert .cellml to .so unless called collectively, due to possible race conditions.");
        TS_ASSERT_THROWS_CONTAINS(converter.Convert(cellml_file),"Conversion of CellML to Chaste shared object failed.");
        EXPECT0(chdir, "..");
    }
    
    void TestArchiving() throw(Exception)
    {
#ifdef CHASTE_CAN_CHECKPOINT_DLLS
        // Copy CellML file into output dir
        std::string dirname = "TestDynamicallyLoadedCellModelsArchiving";
        OutputFileHandler handler(dirname);
        if (PetscTools::AmMaster())
        {
            FileFinder cellml_file("heart/dynamic/luo_rudy_1991.cellml", RelativeTo::ChasteSourceRoot);
            EXPECT0(system, "cp " + cellml_file.GetAbsolutePath() + " " + handler.GetOutputDirectoryFullPath());
        }
        PetscTools::Barrier("TestArchiving_cp");

        // Convert to .so
        CellMLToSharedLibraryConverter converter;
        FileFinder cellml_file(dirname + "/luo_rudy_1991.cellml", RelativeTo::ChasteTestOutput);
        TS_ASSERT(cellml_file.Exists());
        FileFinder so_file(dirname + "/libluo_rudy_1991.so", RelativeTo::ChasteTestOutput);
        TS_ASSERT(!so_file.Exists());
        DynamicCellModelLoader* p_loader = converter.Convert(cellml_file);
        
        // Load a cell model from the .so
        AbstractCardiacCell* p_cell = CreateLr91CellFromLoader(*p_loader, 0u);
        
        // Archive it
        handler.SetArchiveDirectory();
        std::string archive_filename1 = ArchiveLocationInfo::GetProcessUniqueFilePath("first-save.arch");
        {
            AbstractCardiacCell* const p_const_cell = p_cell;
            std::ofstream ofs(archive_filename1.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_const_cell;
        }
        
        // Load from archive
        AbstractCardiacCell* p_loaded_cell1;
        {
            std::ifstream ifs(archive_filename1.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_loaded_cell1;
        }
        
        // Archive the un-archived model
        std::string archive_filename2 = ArchiveLocationInfo::GetProcessUniqueFilePath("second-save.arch");
        {
            AbstractCardiacCell* const p_const_cell = p_loaded_cell1;
            std::ofstream ofs(archive_filename2.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_const_cell;
        }
            
        // Load from the new archive
        AbstractCardiacCell* p_loaded_cell2;
        {
            std::ifstream ifs(archive_filename2.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_loaded_cell2;
        }
        
        // Check simulations of both loaded cells
        SimulateLr91AndCompare(p_loaded_cell1);
        delete p_loaded_cell1;
        SimulateLr91AndCompare(p_loaded_cell2);
        delete p_loaded_cell2;
        delete p_cell;
#else
        std::cout << "Note: this test can only actually test anything on Boost>=1.37." << std::endl;
#endif // CHASTE_CAN_CHECKPOINT_DLLS
    }
};

#endif /* TESTDYNAMICALLYLOADEDCELLMODELS_HPP_ */
