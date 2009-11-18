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


#ifndef TESTHDF5TOVISUALIZERCONVERTERS_HPP_
#define TESTHDF5TOVISUALIZERCONVERTERS_HPP_

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "Hdf5ToMeshalyzerConverter.hpp"
#include "Hdf5ToCmguiConverter.hpp"
#include "Hdf5ToVtkConverter.hpp"
#include "PetscTools.hpp"
#include "OutputFileHandler.hpp"
#include "HeartConfig.hpp"
#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"

typedef Hdf5ToVtkConverter<3,3> VTK_3D;
typedef Hdf5ToCmguiConverter<3,3> CMGUI_3D;
typedef Hdf5ToMeshalyzerConverter<3,3> MESHA_3D;

class TestHdf5ToVisualizerConverters : public CxxTest::TestSuite
{
private :
    // copies a file (relative to Chaste home to CHASTE_TEST_OUTPUT/dir
    void CopyToTestOutputDirectory(std::string file, std::string dir)
    {
        if (PetscTools::AmMaster())
        {
            std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
            std::string command = "mkdir -p " + test_output_directory + dir;
            int return_value;
            return_value = system(command.c_str());
            assert(return_value==0);
            command = "cp " + file + " " + test_output_directory + dir+"/";
            return_value = system(command.c_str());
            assert(return_value==0);
        }
        PetscTools::Barrier();
    }

public :
    void TestMonodomainMeshalyzerConversion() throw(Exception)
    {
        OutputFileHandler handler("TestHdf5ToMeshalyzerConverter");

        // firstly, copy ./heart/test/data/MonoDg01d/*.h5 to CHASTE_TEST_OUTPUT/TestHdf5ToMeshalyzerConverter,
        // as that is where the reader reads from.
        CopyToTestOutputDirectory("heart/test/data/Monodomain1d/MonodomainLR91_1d.h5",
                                  "TestHdf5ToMeshalyzerConverter");

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_100_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // convert
        HeartConfig::Instance()->SetOutputDirectory("TestHdf5ToMeshalyzerConverter");
        Hdf5ToMeshalyzerConverter<1,1> converter("TestHdf5ToMeshalyzerConverter", "MonodomainLR91_1d", &mesh);

        // compare the voltage file with a correct version
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        std::string command = "cmp " + test_output_directory + "/TestHdf5ToMeshalyzerConverter/output/MonodomainLR91_1d_V.dat "
                                     + "heart/test/data/Monodomain1d/MonodomainLR91_1d_V.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);
        command = "cmp " + test_output_directory + "/TestHdf5ToMeshalyzerConverter/output/MonodomainLR91_1d_times.info "
                                     + "heart/test/data/Monodomain1d/MonodomainLR91_1d_times.info";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);
    }

    void TestBidomainMeshalyzerConversion() throw(Exception)
    {
        OutputFileHandler handler("TestHdf5ToMeshalyzerConverter");

        // firstly, copy ./heart/test/data/Bidomain1d/*.h5 to CHASTE_TEST_OUTPUT/TestHdf5ToMeshalyzerConverter,
        // as that is where the reader reads from.
        CopyToTestOutputDirectory("heart/test/data/Bidomain1d/bidomain.h5",
                                  "TestHdf5ToMeshalyzerConverter");

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_100_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // convert
        HeartConfig::Instance()->SetOutputDirectory("TestHdf5ToMeshalyzerConverter");
        Hdf5ToMeshalyzerConverter<1,1> converter("TestHdf5ToMeshalyzerConverter",  "bidomain", &mesh);

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


    void TestMonodomainCmguiConversion3D() throw(Exception)
    {
        std::string working_directory = "TestHdf5ToCmguiConverter_monodomain";
        OutputFileHandler handler(working_directory);

        // firstly, copy ./heart/test/data/CmguiData/*.h5 to CHASTE_TEST_OUTPUT/TestHdf5ToCmguiConverter_monodomain,
        // as that is where the reader reads from.
        CopyToTestOutputDirectory("heart/test/data/CmguiData/monodomain/cube_2mm_12_elements.h5",
                                  working_directory);

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements"); //Not used in the test for exceptions
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // convert
        HeartConfig::Instance()->SetOutputDirectory(working_directory);
        Hdf5ToCmguiConverter<3,3> converter(working_directory, "cube_2mm_12_elements", &mesh);

        // compare the voltage file with a correct version
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        std::string command_first_time_step = "cmp " + test_output_directory + working_directory +"/cmgui_output/cube_2mm_12_elements_0.exnode"
                                     + " heart/test/data/CmguiData/monodomain/cube_2mm_12_elements_0.exnode";
        TS_ASSERT_EQUALS(system(command_first_time_step.c_str()), 0);

        std::string command_second_time_step = "cmp " + test_output_directory + working_directory +"/cmgui_output/cube_2mm_12_elements_1.exnode"
                                     + " heart/test/data/CmguiData/monodomain/cube_2mm_12_elements_1.exnode";
        TS_ASSERT_EQUALS(system(command_second_time_step.c_str()), 0);
    }


    void TestBidomainCmguiConversion3D() throw(Exception)
    {
        std::string working_directory = "TestHdf5ToCmguiConverter_bidomain";
        OutputFileHandler handler(working_directory);

        // firstly, copy ./heart/test/data/CmguiData/*.h5 to CHASTE_TEST_OUTPUT/TestHdf5ToCmguiConverter_monodomain,
        // as that is where the reader reads from.
        CopyToTestOutputDirectory("heart/test/data/CmguiData/bidomain/cube_2mm_12_elements.h5",
                                  working_directory);

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements"); //Not used in the test for exceptions
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // convert
        HeartConfig::Instance()->SetOutputDirectory(working_directory);
        Hdf5ToCmguiConverter<3,3> converter(working_directory, "cube_2mm_12_elements", &mesh);

        // compare the voltage file with a correct version that is known to visualize correctly in Cmgui
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        std::string command_first_time_step = "cmp " + test_output_directory + working_directory +"/cmgui_output/cube_2mm_12_elements_0.exnode"
                                     + " heart/test/data/CmguiData/bidomain/cube_2mm_12_elements_0.exnode";
        TS_ASSERT_EQUALS(system(command_first_time_step.c_str()), 0);

        std::string command_second_time_step = "cmp " + test_output_directory + working_directory +"/cmgui_output/cube_2mm_12_elements_1.exnode"
                                     + " heart/test/data/CmguiData/bidomain/cube_2mm_12_elements_1.exnode";
        TS_ASSERT_EQUALS(system(command_second_time_step.c_str()), 0);
    }

    void TestMonodomainCmguiConversion2D() throw(Exception)
    {
        std::string working_directory = "TestHdf5ToCmguiConverter_monodomain2D";
        OutputFileHandler handler(working_directory);

        // firstly, copy ./heart/test/data/CmguiData/*.h5 to CHASTE_TEST_OUTPUT/TestHdf5ToCmguiConverter_monodomain2D,
        // as that is where the reader reads from. This data file was generated
        // on this mesh by TestMonodomainProblem2DWithPointStimulusInTheVeryCentreOfTheMesh
        CopyToTestOutputDirectory("heart/test/data/CmguiData/monodomain/2D_0_to_1mm_400_elements.h5",
                                  working_directory);

        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_400_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // convert
        HeartConfig::Instance()->SetOutputDirectory(working_directory);
        Hdf5ToCmguiConverter<2,2> converter(working_directory, "2D_0_to_1mm_400_elements", &mesh);

        // compare the voltage file with a correct version that visualizes Vm correctly in cmgui
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        std::string command_first_time_step = "cmp " + test_output_directory + working_directory +"/cmgui_output/2D_0_to_1mm_400_elements_0.exnode"
                                     + " heart/test/data/CmguiData/monodomain/2D_0_to_1mm_400_elements_0.exnode";
        TS_ASSERT_EQUALS(system(command_first_time_step.c_str()), 0);

        std::string command_second_time_step = "cmp " + test_output_directory + working_directory +"/cmgui_output/2D_0_to_1mm_400_elements_1.exnode"
                                     + " heart/test/data/CmguiData/monodomain/2D_0_to_1mm_400_elements_1.exnode";
        TS_ASSERT_EQUALS(system(command_second_time_step.c_str()), 0);
    }

    void TestBidomainCmguiConversion1D() throw(Exception)
    {
        std::string working_directory = "TestHdf5ToCmguiConverter_bidomain1D";
        OutputFileHandler handler(working_directory);

        // firstly, copy ./heart/test/data/CmguiData/*.h5 to CHASTE_TEST_OUTPUT/TestHdf5ToCmguiConverter_bidomain1D,
        // as that is where the reader reads from.
        CopyToTestOutputDirectory("heart/test/data/CmguiData/bidomain/1D_0_to_1_100_elements.h5",
                                  working_directory);

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_100_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // convert
        HeartConfig::Instance()->SetOutputDirectory(working_directory);
        Hdf5ToCmguiConverter<1,1> converter(working_directory, "1D_0_to_1_100_elements", &mesh);

        // compare the voltage file with a correct version that visualizes both Vm and Phie correctly in cmgui
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        std::string command_first_time_step = "cmp " + test_output_directory + working_directory +"/cmgui_output/1D_0_to_1_100_elements_0.exnode"
                                     + " heart/test/data/CmguiData/bidomain/1D_0_to_1_100_elements_0.exnode";
        TS_ASSERT_EQUALS(system(command_first_time_step.c_str()), 0);

        std::string command_second_time_step = "cmp " + test_output_directory + working_directory +"/cmgui_output/1D_0_to_1_100_elements_1.exnode"
                                     + " heart/test/data/CmguiData/bidomain/1D_0_to_1_100_elements_1.exnode";
        TS_ASSERT_EQUALS(system(command_second_time_step.c_str()), 0);
    }


    void TestBidomainVtkConversion3D() throw(Exception)
    {
#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        std::string working_directory = "TestHdf5ToVtkConverter_bidomain";
        OutputFileHandler handler(working_directory);

        // firstly, copy ./heart/test/data/CmguiData/*.h5 to CHASTE_TEST_OUTPUT/TestHdf5ToVtkConverter_bidomain,
        // as that is where the reader reads from.
        CopyToTestOutputDirectory("heart/test/data/CmguiData/bidomain/cube_2mm_12_elements.h5",
                                  working_directory);

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // convert
        HeartConfig::Instance()->SetOutputDirectory(working_directory);
        Hdf5ToVtkConverter<3,3> converter(working_directory, "cube_2mm_12_elements", &mesh);

        // compare the voltage file with a correct version that is known to visualize correctly in Vtk
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        std::string command_first_time_step = "cmp " + test_output_directory + working_directory +"/vtk_output/cube_2mm_12_elements.vtu"
                                     + " heart/test/data/VtkData/bidomain/cube_2mm_12_elements.vtu";
        TS_ASSERT_EQUALS(system(command_first_time_step.c_str()), 0);
#endif //CHASTE_VTK
    }
    void TestMonodomainVtkConversion2D() throw(Exception)
    {
#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        std::string working_directory = "TestHdf5ToVtkConverter_monodomain2D";
        OutputFileHandler handler(working_directory);

        CopyToTestOutputDirectory("heart/test/data/CmguiData/monodomain/2D_0_to_1mm_400_elements.h5",
                                  working_directory);

        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_400_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // convert
        HeartConfig::Instance()->SetOutputDirectory(working_directory);
        Hdf5ToVtkConverter<2,2> converter(working_directory, "2D_0_to_1mm_400_elements", &mesh);

        // compare the voltage file with a correct version that visualizes Vm correctly in VTK
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        std::string command_first_time_step = "cmp " + test_output_directory + working_directory +"/vtk_output/2D_0_to_1mm_400_elements.vtu"
                                     + " heart/test/data/VtkData/monodomain/2D_0_to_1mm_400_elements.vtu";
        TS_ASSERT_EQUALS(system(command_first_time_step.c_str()), 0);
#endif //CHASTE_VTK

    }

    void TestExceptions() throw(Exception)
    {
        std::string directory = "TestHdf5ConverterExceptions";

        CopyToTestOutputDirectory("io/test/data/hdf5_test_full_format.h5", // doesn't have one or two variables
                                  directory);

        HeartConfig::Instance()->SetOutputDirectory(directory);

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements"); //Not used in the test for exceptions until number of nodes is checked
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_THROWS_THIS( MESHA_3D converter(directory, "hdf5_test_full_format", &mesh),
                "Data has zero or more than two variables - doesn\'t appear to be mono or bidomain");

        CopyToTestOutputDirectory("heart/test/data/bad_heart_data_1.h5", // monodomain, with "Volt" instead of "V"
                                  directory);

        TS_ASSERT_THROWS_THIS( CMGUI_3D converter(directory, "bad_heart_data_1", &mesh),
                "One variable, but it is not called \'V\'");

        CopyToTestOutputDirectory("heart/test/data/bad_heart_data_2.h5", // bidomain, with "Volt" instead of "V"
                                  directory);

        TS_ASSERT_THROWS_THIS( VTK_3D converter(directory, "bad_heart_data_2", &mesh),
                "Two variables, but they are not called \'V\' and \'Phi_e\'");

        CopyToTestOutputDirectory("heart/test/data/bad_heart_data_2.h5", // bidomain, with "Volt" instead of "V"
                                  directory);

        CopyToTestOutputDirectory("heart/test/data/CmguiData/monodomain/2D_0_to_1mm_400_elements.h5",
                                  directory);
        TS_ASSERT_THROWS_THIS( CMGUI_3D converter(directory, "2D_0_to_1mm_400_elements", &mesh),
                "Mesh and HDF5 file have a different number of nodes");
    }

};
#endif /*TESTHDF5TOVISUALIZERCONVERTERS_HPP_*/
