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


#ifndef TESTHDF5TOMCMGUICONVERTER_HPP_
#define TESTHDF5TOCMGUICONVERTER_HPP_

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "Hdf5ToCmguiConverter.hpp"
#include "PetscTools.hpp"
#include "OutputFileHandler.hpp"
#include "HeartConfig.hpp"

class TestHdf5ToCmguiConverter : public CxxTest::TestSuite
{
private :
    // copies a file (relative to Chaste home to CHASTE_TEST_OUTPUT/dir)
    void CopyToTestOutputDirectory(std::string file, std::string dir)
    {
        if (PetscTools::AmMaster())
        {
            std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
            std::string command = "cp " + file + " " + test_output_directory + dir+"/";
            int return_value;
            return_value = system(command.c_str());
            assert(return_value==0);
        }
        MPI_Barrier(PETSC_COMM_WORLD);
    }

public :
    void TestMonodomainConversion() throw(Exception)
    {
        std::string working_directory = "TestHdf5ToCmguiConverter_monodomain";
        OutputFileHandler handler(working_directory);

        // firstly, copy ./heart/test/data/CmguiData/*.h5 to CHASTE_TEST_OUTPUT/TestHdf5ToCmguiConverter_monodomain,
        // as that is where the reader reads from.
        CopyToTestOutputDirectory("heart/test/data/CmguiData/monodomain/cube_2mm_12_elements.h5",
                                  working_directory);

        // convert
        HeartConfig::Instance()->SetOutputDirectory(working_directory);
        Hdf5ToCmguiConverter converter(working_directory, "cube_2mm_12_elements");

        // compare the voltage file with a correct version
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        std::string command_first_time_step = "cmp " + test_output_directory + working_directory +"/output/cube_2mm_12_elements_0.exnode"
                                     + " heart/test/data/CmguiData/monodomain/cube_2mm_12_elements_0.exnode";
        TS_ASSERT_EQUALS(system(command_first_time_step.c_str()), 0);

        std::string command_second_time_step = "cmp " + test_output_directory + working_directory +"/output/cube_2mm_12_elements_1.exnode"
                                     + " heart/test/data/CmguiData/monodomain/cube_2mm_12_elements_1.exnode";
        TS_ASSERT_EQUALS(system(command_second_time_step.c_str()), 0);
    }


    void TestBidomainConversion() throw(Exception)
    {
        std::string working_directory = "TestHdf5ToCmguiConverter_bidomain";
        OutputFileHandler handler(working_directory);

        // firstly, copy ./heart/test/data/CmguiData/*.h5 to CHASTE_TEST_OUTPUT/TestHdf5ToCmguiConverter_monodomain,
        // as that is where the reader reads from.
        CopyToTestOutputDirectory("heart/test/data/CmguiData/bidomain/cube_2mm_12_elements.h5",
                                  working_directory);

        // convert
        HeartConfig::Instance()->SetOutputDirectory(working_directory);
        Hdf5ToCmguiConverter converter(working_directory, "cube_2mm_12_elements");

        // compare the voltage file with a correct version
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        std::string command_first_time_step = "cmp " + test_output_directory + working_directory +"/output/cube_2mm_12_elements_0.exnode"
                                     + " heart/test/data/CmguiData/bidomain/cube_2mm_12_elements_0.exnode";
        TS_ASSERT_EQUALS(system(command_first_time_step.c_str()), 0);

        std::string command_second_time_step = "cmp " + test_output_directory + working_directory +"/output/cube_2mm_12_elements_1.exnode"
                                     + " heart/test/data/CmguiData/bidomain/cube_2mm_12_elements_1.exnode";
        TS_ASSERT_EQUALS(system(command_second_time_step.c_str()), 0);
    }

    void TestExceptions() throw(Exception)
    {
        std::string bidomain_directory = "TestHdf5ToCmguiConverter_bidomain";
        std::string monodomain_directory = "TestHdf5ToCmguiConverter_monodomain";

        CopyToTestOutputDirectory("io/test/data/hdf5_test_full_format.h5", // doesn't have one or two variables
                                  bidomain_directory);

        HeartConfig::Instance()->SetOutputDirectory(bidomain_directory);

        TS_ASSERT_THROWS_THIS( Hdf5ToCmguiConverter converter(bidomain_directory, "hdf5_test_full_format"),
                "Data has zero or more than two variables - doesn\'t appear to be mono or bidomain");

        CopyToTestOutputDirectory("heart/test/data/bad_heart_data_1.h5", // monodomain, with "Volt" instead of "V"
                                  monodomain_directory);

        TS_ASSERT_THROWS_THIS( Hdf5ToCmguiConverter converter2("TestHdf5ToCmguiConverter_monodomain", "bad_heart_data_1"),
                "One variable, but it is not called \'V\'");

        CopyToTestOutputDirectory("heart/test/data/bad_heart_data_2.h5", // bidomain, with "Volt" instead of "V"
                                  bidomain_directory);

        TS_ASSERT_THROWS_THIS( Hdf5ToCmguiConverter converter2(bidomain_directory, "bad_heart_data_2"),
                "Two variables, but they are not called \'V\' and \'Phi_e\'");
    }
};
#endif /*TESTHDF5TOCMGUICONVERTER_HPP_*/
