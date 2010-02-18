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
#include "HeartConfig.hpp"
#include "FileFinder.hpp"
#include "CellMLToSharedLibraryConverter.hpp"

class TestDynamicallyLoadedCellModels : public CxxTest::TestSuite
{
private:

    void RunLr91Test(DynamicCellModelLoader& rLoader,
                     unsigned v_index=4u)
    {
        // Set stimulus
        double magnitude = -25.5;
        double duration  = 2.0 ;  // ms
        double when = 50.0; // ms

        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude, duration, when));
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        double end_time = 1000.0; //One second in milliseconds
        
        // Load the cell model dynamically
        AbstractCardiacCell* p_cell = rLoader.CreateCell(p_solver, p_stimulus);

        // Simple sanity check
        TS_ASSERT_EQUALS(p_cell->GetVoltageIndex(), v_index);

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
public:
    /**
     * This is based on TestOdeSolverForLR91WithDelayedSimpleStimulus from
     * TestIonicModels.hpp.
     */
    void TestDynamicallyLoadedLr91() throw(Exception)
    {
        // Load the cell model dynamically
        std::string model_name = "libDynamicallyLoadableLr91.so";
        DynamicCellModelLoader loader(ChasteComponentBuildDir("heart") + "dynamic/" + model_name);
        RunLr91Test(loader);
        
        // The .so also gets copied into the source folder
        DynamicCellModelLoader loader2(std::string(ChasteBuildRootDir()) + "heart/dynamic/" + model_name);
        RunLr91Test(loader2);
    }
    
    /**
     * Currently the build system will not automatically regenerate the C++ code from CellML, so you
     * need to do that step manually:
     *     cdchaste
     *     ./python/ConvertCellModel.py -y heart/dynamic/luo_rudy_1991.cellml
     */
    void TestLr91FromCellML() throw(Exception)
    {
        FileFinder model("heart/dynamic/libluo_rudy_1991.so", cp::relative_to_type::chaste_source_root);
        DynamicCellModelLoader loader(model.GetAbsolutePath());
        RunLr91Test(loader, 0u);
        
        FileFinder model_opt("heart/dynamic/libluo_rudy_1991Opt.so", cp::relative_to_type::chaste_source_root);
        DynamicCellModelLoader loader_opt(model_opt.GetAbsolutePath());
        RunLr91Test(loader_opt, 0u);
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
        
        // Now mock up what HeartConfigRelatedCellFactory will have to do
        TS_ASSERT(ionic_model.Dynamic().present());
        if (ionic_model.Dynamic().present())
        {
            FileFinder file_finder(ionic_model.Dynamic()->Path());
            TS_ASSERT(file_finder.Exists());
            DynamicCellModelLoader loader(file_finder.GetAbsolutePath());
            
            RunLr91Test(loader);
        }
        // then 'else' what it currently does, more or less...
    }
    
    void TestCellmlConverter() throw(Exception)
    {
        // Copy CellML file into output dir
        std::string dirname = "TestCellmlConverter";
        OutputFileHandler handler(dirname);
        if (PetscTools::AmMaster())
        {
            FileFinder cellml_file("heart/dynamic/luo_rudy_1991.cellml", cp::relative_to_type::chaste_source_root);
            EXPECT0(system, "cp " + cellml_file.GetAbsolutePath() + " " + handler.GetOutputDirectoryFullPath());
        }
        PetscTools::Barrier("TestCellmlConverter");
        
        CellMLToSharedLibraryConverter converter;
        
        // Convert a real CellML file
        FileFinder cellml_file(dirname + "/luo_rudy_1991.cellml", cp::relative_to_type::chaste_test_output);
        TS_ASSERT(cellml_file.Exists());
        FileFinder so_file(dirname + "/libluo_rudy_1991.so", cp::relative_to_type::chaste_test_output);
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
        FileFinder no_ext(dirname + "/" + file_name, cp::relative_to_type::chaste_test_output);
        TS_ASSERT_THROWS_THIS(converter.Convert(no_ext), "Dynamically loadable cell model '"
                              + no_ext.GetAbsolutePath() + "' does not exist.");
        
        out_stream fp = handler.OpenOutputFile(file_name);
        fp->close();
        TS_ASSERT_THROWS_THIS(converter.Convert(no_ext), "File does not have an extension: " + no_ext.GetAbsolutePath());
        FileFinder unsupp_ext("heart/src/io/FileFinder.hpp", cp::relative_to_type::chaste_source_root);
        TS_ASSERT_THROWS_THIS(converter.Convert(unsupp_ext), "Unsupported extension '.hpp' of file '"
                              + unsupp_ext.GetAbsolutePath() + "'; must be .so or .cellml");

        EXPECT0(chdir, "heart"); // The ConvertCellModel.py script in ConvertCellmlToSo() should only work from chaste source directory.
        // Having the 'rm' after the 'chdir' cunningly double-checks that the FileFinder has really given us an absolute path
        EXPECT0(system, "rm " + so_file.GetAbsolutePath()); // Make sure the conversion is re-run
        TS_ASSERT_THROWS_CONTAINS(converter.Convert(cellml_file),"Conversion of CellML to Chaste shared object failed.");
        EXPECT0(chdir, "..");
    }
};

#endif /* TESTDYNAMICALLYLOADEDCELLMODELS_HPP_ */
