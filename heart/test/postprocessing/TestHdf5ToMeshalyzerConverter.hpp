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


#ifndef TESTHDF5TOMESHALYZERCONVERTER_HPP_
#define TESTHDF5TOMESHALYZERCONVERTER_HPP_

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "Hdf5ToMeshalyzerConverter.hpp"
#include "PetscTools.hpp"
#include "OutputFileHandler.hpp"
#include "HeartConfig.hpp"

class TestHdf5ToMeshalyzerConverter : public CxxTest::TestSuite
{
private :
    // copies a file (relative to Chaste home to CHASTE_TEST_OUTPUT/dir
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
        OutputFileHandler handler("TestHdf5ToMeshalyzerConverter");

        // firstly, copy ./heart/test/data/MonoDg01d/*.h5 to CHASTE_TEST_OUTPUT/TestHdf5ToMeshalyzerConverter,
        // as that is where the reader reads from.
        CopyToTestOutputDirectory("heart/test/data/Monodomain1d/MonodomainLR91_1d.h5",
                                  "TestHdf5ToMeshalyzerConverter");

        // convert
        HeartConfig::Instance()->SetOutputDirectory("TestHdf5ToMeshalyzerConverter");
        Hdf5ToMeshalyzerConverter converter("TestHdf5ToMeshalyzerConverter", "MonodomainLR91_1d");

        // compare the voltage file with a correct version
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        std::string command = "cmp " + test_output_directory + "/TestHdf5ToMeshalyzerConverter/output/MonodomainLR91_1d_V.dat "
                                     + "heart/test/data/Monodomain1d/MonodomainLR91_1d_V.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);
        command = "cmp " + test_output_directory + "/TestHdf5ToMeshalyzerConverter/output/MonodomainLR91_1d_times.info "
                                     + "heart/test/data/Monodomain1d/MonodomainLR91_1d_times.info";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);
        }


    void TestBidomainConversion() throw(Exception)
    {
        OutputFileHandler handler("TestHdf5ToMeshalyzerConverter");

        // firstly, copy ./heart/test/data/Bidomain1d/*.h5 to CHASTE_TEST_OUTPUT/TestHdf5ToMeshalyzerConverter,
        // as that is where the reader reads from.
        CopyToTestOutputDirectory("heart/test/data/Bidomain1d/bidomain.h5",
                                  "TestHdf5ToMeshalyzerConverter");

        // convert
        HeartConfig::Instance()->SetOutputDirectory("TestHdf5ToMeshalyzerConverter");
        Hdf5ToMeshalyzerConverter converter("TestHdf5ToMeshalyzerConverter",  "bidomain");

        // compare the voltage file
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        std::string command = "cmp " + test_output_directory + "/TestHdf5ToMeshalyzerConverter/output/bidomain_V.dat "
                                     + "heart/test/data/Bidomain1d/bidomain_V.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);

        // compare the Phi_e file
        command = "cmp " + test_output_directory + "/TestHdf5ToMeshalyzerConverter/output/bidomain_Phi_e.dat "
                         + "heart/test/data/Bidomain1d/bidomain_Phi_e.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);

       // compare the time information file
        command = "cmp " + test_output_directory + "/TestHdf5ToMeshalyzerConverter/output/bidomain_times.info "
                         + "heart/test/data/Bidomain1d/bidomain_times.info";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);
    }

    void TestExceptions() throw(Exception)
    {
        OutputFileHandler handler("TestHdf5ToMeshalyzerConverter");

        CopyToTestOutputDirectory("io/test/data/hdf5_test_full_format.h5", // doesn't have one or two variables
                                  "TestHdf5ToMeshalyzerConverter");

        HeartConfig::Instance()->SetOutputDirectory("TestHdf5ToMeshalyzerConverter");

        TS_ASSERT_THROWS_THIS( Hdf5ToMeshalyzerConverter converter("TestHdf5ToMeshalyzerConverter", "hdf5_test_full_format"),
                "Data has zero or more than two variables - doesn\'t appear to be mono or bidomain");

        CopyToTestOutputDirectory("heart/test/data/bad_heart_data_1.h5", // monodomain, with "Volt" instead of "V"
                                  "TestHdf5ToMeshalyzerConverter");

        TS_ASSERT_THROWS_THIS( Hdf5ToMeshalyzerConverter converter2("TestHdf5ToMeshalyzerConverter", "bad_heart_data_1"),
                "One variable, but it is not called \'V\'");

        CopyToTestOutputDirectory("heart/test/data/bad_heart_data_2.h5", // bidomain, with "Volt" instead of "V"
                                  "TestHdf5ToMeshalyzerConverter");

        TS_ASSERT_THROWS_THIS( Hdf5ToMeshalyzerConverter converter2("TestHdf5ToMeshalyzerConverter", "bad_heart_data_2"),
                "Two variables, but they are not called \'V\' and \'Phi_e\'");
    }
};
#endif /*TESTHDF5TOMESHALYZERCONVERTER_HPP_*/
