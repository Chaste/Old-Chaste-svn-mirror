/*
Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef TESTHDF5TOMESHALYZERCONVERTER_HPP_
#define TESTHDF5TOMESHALYZERCONVERTER_HPP_

#include <cxxtest/TestSuite.h>
#include "PetscSetupAndFinalize.hpp"
#include "Hdf5ToMeshalyzerConverter.hpp"

class TestHdf5ToMeshalyzerConverter : public CxxTest::TestSuite
{
private :
    void CopyToTestOutputDirectory(std::string file, std::string dir)
    {
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        std::string command = "cp " + file + " " + test_output_directory + dir;
        int return_value = system(command.c_str());
        assert(return_value==0);
    }

public :
    
    void TestMonodomainConvertion() throw(Exception)
    {
        OutputFileHandler handler("TestHdf5ToMeshalyzerConverter");
        
        // firstly, copy ./heart/test/data/MonoDg01d/*.h5 to CHASTE_TEST_OUTPUT/TestHdf5ToMeshalyzerConverter,
        // as that is where the reader reads from.
        CopyToTestOutputDirectory("heart/test/data/Monodomain1d/MonodomainLR91_1d.h5",
                                  "TestHdf5ToMeshalyzerConverter");

        // convert
        Hdf5ToMeshalyzerConverter converter("TestHdf5ToMeshalyzerConverter", "MonodomainLR91_1d");
        
        // compare the voltage file
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        std::string command = "cmp " + test_output_directory + "/TestHdf5ToMeshalyzerConverter/MonodomainLR91_1d_V.dat " 
                                     + "heart/test/data/Monodomain1d/MonodomainLR91_1d_V.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);
    }


    void TestBidomainConvertion() throw(Exception)
    {
        OutputFileHandler handler("TestHdf5ToMeshalyzerConverter");

        // firstly, copy ./heart/test/data/Bidomain1d/*.h5 to CHASTE_TEST_OUTPUT/TestHdf5ToMeshalyzerConverter,
        // as that is where the reader reads from.
        CopyToTestOutputDirectory("heart/test/data/Bidomain1d/bidomain.h5",
                                  "TestHdf5ToMeshalyzerConverter");

        // convert
        Hdf5ToMeshalyzerConverter converter("TestHdf5ToMeshalyzerConverter", "bidomain");

        // compare the voltage file
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        std::string command = "cmp " + test_output_directory + "/TestHdf5ToMeshalyzerConverter/bidomain_V.dat " 
                                     + "heart/test/data/Bidomain1d/bidomain_V.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);

        // compare the Phi_e file
        command = "cmp " + test_output_directory + "/TestHdf5ToMeshalyzerConverter/bidomain_Phi_e.dat " 
                         + "heart/test/data/Bidomain1d/bidomain_Phi_e.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);
    }
    
    void TestExceptions() throw(Exception)
    {
        OutputFileHandler handler("TestHdf5ToMeshalyzerConverter");

        CopyToTestOutputDirectory("io/test/data/hdf5_test_full_format.h5", // doesn't have one two variables
                                  "TestHdf5ToMeshalyzerConverter");

        TS_ASSERT_THROWS_ANYTHING( Hdf5ToMeshalyzerConverter converter("TestHdf5ToMeshalyzerConverter", "hdf5_test_full_format") );

        CopyToTestOutputDirectory("heart/test/data/bad_heart_data_1.h5", // monodomain, with "Volt" instead of "V"
                                  "TestHdf5ToMeshalyzerConverter");

        TS_ASSERT_THROWS_ANYTHING( Hdf5ToMeshalyzerConverter converter2("TestHdf5ToMeshalyzerConverter", "bad_heart_data_1") );


        CopyToTestOutputDirectory("heart/test/data/bad_heart_data_2.h5", // bidomain, with "Volt" instead of "V"
                                  "TestHdf5ToMeshalyzerConverter");

        TS_ASSERT_THROWS_ANYTHING( Hdf5ToMeshalyzerConverter converter2("TestHdf5ToMeshalyzerConverter", "bad_heart_data_2") );

    }
};
#endif /*TESTHDF5TOMESHALYZERCONVERTER_HPP_*/
