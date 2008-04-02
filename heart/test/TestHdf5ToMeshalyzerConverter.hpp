#ifndef TESTHDF5TOMESHALYZERCONVERTER_HPP_
#define TESTHDF5TOMESHALYZERCONVERTER_HPP_

#include <cxxtest/TestSuite.h>
#include "PetscSetupAndFinalize.hpp"
#include "Hdf5ToMeshalyzerConverter.hpp"

class TestHdf5ToMeshalyzerConverter : public CxxTest::TestSuite
{
private:
    // the setup method just creates the test directory that will be used 
    void setUpWorld()
    {
        OutputFileHandler handler("TestHdf5ToMeshalyzerConverter");
    }

public:
    
    void TestMonodomainConvertion() throw(Exception)
    {
        // firstly, copy ./heart/test/data/MonoDg01d/*.h5 to CHASTE_TEST_OUTPUT/TestHdf5ToMeshalyzerConverter,
        // as that is where the reader reads from.
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        std::string command = "cp heart/test/data/MonoDg01d/NewMonodomainLR91_1d.h5 " + test_output_directory +  "/TestHdf5ToMeshalyzerConverter";
        int return_value = system(command.c_str());
        assert(return_value==0);
        
        // convert
        std::string h5_output_dir = "TestHdf5ToMeshalyzerConverter";
        std::string h5_file_name = "NewMonodomainLR91_1d";
        
        Hdf5ToMeshalyzerConverter converter(h5_output_dir, h5_file_name);
        
        // compare the voltage file
        command = "cmp " + test_output_directory + "/TestHdf5ToMeshalyzerConverter/NewMonodomainLR91_1d_V.dat " 
                         + "heart/test/data/MonoDg01d/NewMonodomainLR91_1d_V.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);
    }


    void TestBidomainConvertion() throw(Exception)
    {
        // firstly, copy ./heart/test/data/Bidomain1d/*.h5 to CHASTE_TEST_OUTPUT/TestHdf5ToMeshalyzerConverter,
        // as that is where the reader reads from.
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        std::string command = "cp heart/test/data/Bidomain1d/bidomain.h5 " + test_output_directory +  "/TestHdf5ToMeshalyzerConverter";
        int return_value = system(command.c_str());
        assert(return_value==0);

        // convert
        std::string h5_output_dir = "TestHdf5ToMeshalyzerConverter";
        std::string h5_file_name = "bidomain";
        
        Hdf5ToMeshalyzerConverter converter(h5_output_dir, h5_file_name);

        // compare the voltage file
        command = "cmp " + test_output_directory + "/TestHdf5ToMeshalyzerConverter/bidomain_V.dat " 
                         + "heart/test/data/Bidomain1d/bidomain_V.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);

        // compare the Phi_e file
        command = "cmp " + test_output_directory + "/TestHdf5ToMeshalyzerConverter/bidomain_Phi_e.dat " 
                         + "heart/test/data/Bidomain1d/bidomain_Phi_e.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);

    }
};
#endif /*TESTHDF5TOMESHALYZERCONVERTER_HPP_*/
