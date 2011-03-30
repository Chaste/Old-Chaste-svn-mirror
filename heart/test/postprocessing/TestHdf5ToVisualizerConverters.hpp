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
#include "DistributedTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"

#ifdef CHASTE_VTK
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkVersion.h>
#endif

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
        std::string command = "diff -a -I \"Created by Chaste\" " + test_output_directory + "/TestHdf5ToMeshalyzerConverter/output/MonodomainLR91_1d_V.dat "
                                     + "heart/test/data/Monodomain1d/MonodomainLR91_1d_V.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);
        command = "diff -a -I \"Created by Chaste\" " + test_output_directory + "/TestHdf5ToMeshalyzerConverter/output/MonodomainLR91_1d_times.info "
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
        std::string command = "diff -a -I \"Created by Chaste\" " + test_output_directory + "/TestHdf5ToMeshalyzerConverter/output/bidomain_V.dat "
                                     + "heart/test/data/Bidomain1d/bidomain_V.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);

        // compare the Phi_e file
        command = "diff -a -I \"Created by Chaste\" " + test_output_directory + "/TestHdf5ToMeshalyzerConverter/output/bidomain_Phi_e.dat "
                         + "heart/test/data/Bidomain1d/bidomain_Phi_e.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);

       // compare the time information file
        command = "diff -a -I \"Created by Chaste\" " + test_output_directory + "/TestHdf5ToMeshalyzerConverter/output/bidomain_times.info "
                         + "heart/test/data/Bidomain1d/bidomain_times.info";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);
    }

    //This test covers the case when the hdf5 file contains 3 variables (e.g., after solving a problem with PROBLEM_DIM=3)
    void TestMeshalyzerConversion3Variables() throw(Exception)
    {
        OutputFileHandler handler("TestMeshalyzerConversion3Variables");

        // firstly, copy ./heart/test/data/Bidomain1d/*.h5 to CHASTE_TEST_OUTPUT/TestHdf5ToMeshalyzerConverter,
        // as that is where the reader reads from.
        CopyToTestOutputDirectory("heart/test/data/three_variables/3_vars.h5",
                                  "TestMeshalyzerConversion3Variables");

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_100_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // convert
        HeartConfig::Instance()->SetOutputDirectory("TestMeshalyzerConversion3Variables");
        Hdf5ToMeshalyzerConverter<1,1> converter("TestMeshalyzerConversion3Variables",  "3_vars", &mesh);

        // compare the first voltage file
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        std::string command = "diff -a -I \"Created by Chaste\" " + test_output_directory + "/TestMeshalyzerConversion3Variables/output/3_vars_Vm_1.dat "
                                     + "heart/test/data/three_variables/extended_bidomain_Vm_1.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);

        // compare the second voltage file

        command = "diff -a -I \"Created by Chaste\" " + test_output_directory + "/TestMeshalyzerConversion3Variables/output/3_vars_Vm_2.dat "
                                     + "heart/test/data/three_variables/extended_bidomain_Vm_2.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);

        // compare the Phi_e file
        command = "diff -a -I \"Created by Chaste\" " + test_output_directory + "/TestMeshalyzerConversion3Variables/output/3_vars_Phi_e.dat "
                         + "heart/test/data/three_variables/extended_bidomain_Phi_e.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);

       // compare the time information file
        command = "diff -a -I \"Created by Chaste\" " + test_output_directory + "/TestMeshalyzerConversion3Variables/output/3_vars_times.info "
                         + "heart/test/data/three_variables/extended_bidomain_times.info";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);
    }
    
    //This test covers the case when the hdf5 file contains more than 3 variables
    void TestMeshalyzerConversionLotsOfVariables() throw(Exception)
    {
        std::string output_dir = "TestHdf5ToMeshalyzerConversionManyVariables";
        OutputFileHandler handler(output_dir);

        // firstly, copy ./heart/test/data/many_variables/*.h5 to CHASTE_TEST_OUTPUT/TestHdf5ToMeshalyzerConverter,
        // as that is where the reader reads from.
        CopyToTestOutputDirectory("heart/test/data/many_variables/many_variables.h5",
                                  output_dir);

        TrianglesMeshReader<1,1> mesh_reader("heart/test/data/many_variables/1D_65_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // convert
        HeartConfig::Instance()->SetOutputDirectory(output_dir);
        Hdf5ToMeshalyzerConverter<1,1> converter(output_dir,  "many_variables", &mesh);

        std::vector<std::string> variable_names;
        variable_names.push_back("V");
        variable_names.push_back("I_ks");
        variable_names.push_back("I_kr");
        variable_names.push_back("I_Ca_tot");
        variable_names.push_back("I_tot");
        variable_names.push_back("I_Na_tot");
               
        
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        std::string command;        
        for (unsigned i=0; i<variable_names.size(); i++)
        {
            // compare the results files
            command = "diff -a -I \"Created by Chaste\" " + test_output_directory + "/" + output_dir + "/output/many_variables_" + variable_names[i] + ".dat "
                                         + "heart/test/data/many_variables/many_variables_" + variable_names[i] + ".dat";
            TS_ASSERT_EQUALS(system(command.c_str()), 0);
        }
       
        // compare the time information file
        command = "diff -a -I \"Created by Chaste\" " + test_output_directory + output_dir + "/output/many_variables_times.info "
                         + "heart/test/data/many_variables/many_variables_times.info";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);
    }
    
    //This test covers the case when the hdf5 file contains more than 3 variables
    void TestCmguiConversionLotsOfVariables() throw(Exception)
    {
        std::string output_dir = "TestHdf5ToCmguiConversionManyVariables";
        OutputFileHandler handler(output_dir);

        // firstly, copy ./heart/test/data/many_variables/*.h5 to CHASTE_TEST_OUTPUT/TestHdf5ToCmguiConverter,
        // as that is where the reader reads from.
        CopyToTestOutputDirectory("heart/test/data/many_variables/many_variables.h5",
                                  output_dir);

        TrianglesMeshReader<1,1> mesh_reader("heart/test/data/many_variables/1D_65_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // convert
        HeartConfig::Instance()->SetOutputDirectory(output_dir);
        Hdf5ToCmguiConverter<1,1> converter(output_dir,  "many_variables", &mesh);
      
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
     
        // compare the results files
        std::string command = "diff -a -I \"Created by Chaste\" " + test_output_directory + "/" + output_dir + "/cmgui_output/many_variables_0.exnode "
                                     + "heart/test/data/many_variables/many_variables_0.exnode";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);

        //check validity of cmgui script
        std::string command_script = "diff -a -I \"Created by Chaste\" " + test_output_directory + output_dir +"/cmgui_output/script.com"
                                     + " heart/test/data/many_variables/CmguiValidScript.com";
        TS_ASSERT_EQUALS(system(command_script.c_str()), 0);
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
        std::string command_first_time_step = "diff -a -I \"Created by Chaste\" " + test_output_directory + working_directory +"/cmgui_output/cube_2mm_12_elements_0.exnode"
                                     + " heart/test/data/CmguiData/monodomain/cube_2mm_12_elements_0.exnode";
        TS_ASSERT_EQUALS(system(command_first_time_step.c_str()), 0);

        std::string command_second_time_step = "diff -a -I \"Created by Chaste\" " + test_output_directory + working_directory +"/cmgui_output/cube_2mm_12_elements_1.exnode"
                                     + " heart/test/data/CmguiData/monodomain/cube_2mm_12_elements_1.exnode";
        TS_ASSERT_EQUALS(system(command_second_time_step.c_str()), 0);

        //check validity of cmgui script
        std::string command_script = "diff -a -I \"Created by Chaste\" " + test_output_directory + working_directory +"/cmgui_output/script.com"
                                     + " heart/test/data/CmguiData/monodomain/monodomain3dValidScript.com";
        TS_ASSERT_EQUALS(system(command_script.c_str()), 0);
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
        std::string command_first_time_step = "diff -a -I \"Created by Chaste\" " + test_output_directory + working_directory +"/cmgui_output/cube_2mm_12_elements_0.exnode"
                                     + " heart/test/data/CmguiData/bidomain/cube_2mm_12_elements_0.exnode";
        TS_ASSERT_EQUALS(system(command_first_time_step.c_str()), 0);

        std::string command_second_time_step = "diff -a -I \"Created by Chaste\" " + test_output_directory + working_directory +"/cmgui_output/cube_2mm_12_elements_1.exnode"
                                     + " heart/test/data/CmguiData/bidomain/cube_2mm_12_elements_1.exnode";
        TS_ASSERT_EQUALS(system(command_second_time_step.c_str()), 0);
    }
    
    void TestBidomainWithBathCmguiConversion1D() throw(Exception)
    {
        std::string working_directory = "TestBidomainWithBathCmguiConversion1D";
        OutputFileHandler handler(working_directory);

        // firstly, copy ./heart/test/data/CmguiData/*.h5 to CHASTE_TEST_OUTPUT/TestHdf5ToCmguiConverter_monodomain,
        // as that is where the reader reads from.
        CopyToTestOutputDirectory("heart/test/data/CmguiData/bidomain_with_bath/bidomain_with_bath_1d.h5",
                                  working_directory);

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements_with_two_attributes"); //Not used in the test for exceptions
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // convert
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain_with_bath_1d");
        HeartConfig::Instance()->SetOutputDirectory(working_directory);
        Hdf5ToCmguiConverter<1,1> converter(working_directory, "bidomain_with_bath_1d", &mesh, true);

        // compare the voltage file with a correct version that is known to visualize correctly in Cmgui
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        
        //mesh file first, one exnode, one for bath and one for tissue
        std::string command_node_file = "diff -a -I \"Created by Chaste\" " + test_output_directory + working_directory +"/cmgui_output/bidomain_with_bath_1d.exnode"
                                     + " heart/test/data/CmguiData/bidomain_with_bath/bidomain_with_bath_1d.exnode";
        TS_ASSERT_EQUALS(system(command_node_file.c_str()), 0);
        
        std::string command_tissue_element_file = "diff -a -I \"Created by Chaste\" " + test_output_directory + working_directory +"/cmgui_output/tissue.exelem"
                                     + " heart/test/data/CmguiData/bidomain_with_bath/tissue.exelem";
        TS_ASSERT_EQUALS(system(command_tissue_element_file.c_str()), 0);

        std::string command_bath_element_file = "diff -a -I \"Created by Chaste\" " + test_output_directory + working_directory +"/cmgui_output/bath.exelem"
                                     + " heart/test/data/CmguiData/bidomain_with_bath/bath.exelem";
        TS_ASSERT_EQUALS(system(command_bath_element_file.c_str()), 0);
        
        //then the data file
        std::string command_first_time_step = "diff -a -I \"Created by Chaste\" " + test_output_directory + working_directory +"/cmgui_output/bidomain_with_bath_1d_0.exnode"
                                     + " heart/test/data/CmguiData/bidomain_with_bath/bidomain_with_bath_1d_0.exnode";
        TS_ASSERT_EQUALS(system(command_first_time_step.c_str()), 0);

        std::string command_second_time_step = "diff -a -I \"Created by Chaste\" " + test_output_directory + working_directory +"/cmgui_output/bidomain_with_bath_1d_1.exnode"
                                     + " heart/test/data/CmguiData/bidomain_with_bath/bidomain_with_bath_1d_1.exnode";
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
        std::string command_first_time_step = "diff -a -I \"Created by Chaste\"  " + test_output_directory + working_directory +"/cmgui_output/2D_0_to_1mm_400_elements_0.exnode"
                                     + " heart/test/data/CmguiData/monodomain/2D_0_to_1mm_400_elements_0.exnode";
        TS_ASSERT_EQUALS(system(command_first_time_step.c_str()), 0);

        std::string command_second_time_step = "diff -a -I \"Created by Chaste\" " + test_output_directory + working_directory +"/cmgui_output/2D_0_to_1mm_400_elements_1.exnode"
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
        std::string command_first_time_step = "diff -a -I \"Created by Chaste\" " + test_output_directory + working_directory +"/cmgui_output/1D_0_to_1_100_elements_0.exnode"
                                     + " heart/test/data/CmguiData/bidomain/1D_0_to_1_100_elements_0.exnode";
        TS_ASSERT_EQUALS(system(command_first_time_step.c_str()), 0);

        std::string command_second_time_step = "diff -a -I \"Created by Chaste\" " + test_output_directory + working_directory +"/cmgui_output/1D_0_to_1_100_elements_1.exnode"
                                     + " heart/test/data/CmguiData/bidomain/1D_0_to_1_100_elements_1.exnode";
        TS_ASSERT_EQUALS(system(command_second_time_step.c_str()), 0);
    }

    void TestCmguiConversion1DWith3Variables() throw(Exception)
    {
        std::string working_directory = "TestHdf5ToCmguiConverter3Vars";
        OutputFileHandler handler(working_directory);

        CopyToTestOutputDirectory("heart/test/data/three_variables/3_vars.h5",
                                  working_directory);

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_100_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // convert
        HeartConfig::Instance()->SetOutputFilenamePrefix("3_vars");
        HeartConfig::Instance()->SetOutputDirectory(working_directory);
        Hdf5ToCmguiConverter<1,1> converter(working_directory, "3_vars", &mesh);

        // compare the voltage file with a correct version that visualizes both Vs and Phie correctly in cmgui
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();

        std::string command_node_mesh_file = "diff -a -I \"Created by Chaste\" " + test_output_directory + working_directory +"/cmgui_output/3_vars.exnode"
                                     + " heart/test/data/CmguiData/extended_bidomain/3_vars.exnode";
        TS_ASSERT_EQUALS(system(command_node_mesh_file.c_str()), 0);

        std::string command_element_mesh_file = "diff -a -I \"Created by Chaste\" " + test_output_directory + working_directory +"/cmgui_output/3_vars.exelem"
                                     + " heart/test/data/CmguiData/extended_bidomain/3_vars.exelem";
        TS_ASSERT_EQUALS(system(command_node_mesh_file.c_str()), 0);

        std::string command_data_file= "diff -a -I \"Created by Chaste\" " + test_output_directory + working_directory +"/cmgui_output/3_vars_25.exnode"
                                     + " heart/test/data/CmguiData/extended_bidomain/3_vars_25.exnode";
        TS_ASSERT_EQUALS(system(command_data_file.c_str()), 0);
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
        Hdf5ToVtkConverter<3,3> converter(working_directory, "cube_2mm_12_elements", &mesh, false);

        //Note that VTK is not thread-safe.  The master process has spawned a child to write the mesh and may still be writing!
        //This barrier just slows things down a bit
        PetscTools::Barrier();

        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        
        VtkMeshReader<3,3> vtk_mesh_reader(test_output_directory + working_directory +"/vtk_output/cube_2mm_12_elements.vtu");
        TS_ASSERT_EQUALS( vtk_mesh_reader.GetNumNodes(), 12U);
        TS_ASSERT_EQUALS( vtk_mesh_reader.GetNumElements(), 12U);
        
        std::vector<double> first_node = vtk_mesh_reader.GetNextNode();
        TS_ASSERT_DELTA( first_node[0] , 0.0 , 1e-6 );
        TS_ASSERT_DELTA( first_node[1] , 0.0, 1e-6 );
        TS_ASSERT_DELTA( first_node[2] , 0.0 , 1e-6 );

        std::vector<double> next_node = vtk_mesh_reader.GetNextNode();
        TS_ASSERT_DELTA( next_node[0] , 0.2 , 1e-6 );
        TS_ASSERT_DELTA( next_node[1] , 0.0 , 1e-6 );
        TS_ASSERT_DELTA( next_node[2] , 0.0 , 1e-6 );
        
        //V_m and phi_e samples 
        std::vector<double> v_at_last, phi_at_last;
        vtk_mesh_reader.GetPointData( "V_000001", v_at_last);
        TS_ASSERT_DELTA( v_at_last[0],  -46.3761, 1e-3 );
        TS_ASSERT_DELTA( v_at_last[6],  -46.3761, 1e-3 );
        TS_ASSERT_DELTA( v_at_last[11], -46.3760, 1e-3 );
        vtk_mesh_reader.GetPointData( "Phi_e_000001", phi_at_last);
        TS_ASSERT_DELTA( phi_at_last[0],  0.0, 1e-3 );
        TS_ASSERT_DELTA( phi_at_last[6],  0.0, 1e-3 );
        TS_ASSERT_DELTA( phi_at_last[11], 0.0, 1e-3 );

        //Show that trying to write pvtu files from a TetrahedralMesh gives a warning (but writes anyway)
        Hdf5ToVtkConverter<3,3> converter2(working_directory, "cube_2mm_12_elements", &mesh, true);
        VtkMeshReader<3,3> vtk_mesh_reader2(test_output_directory + working_directory +"/vtk_output/cube_2mm_12_elements.vtu");
        TS_ASSERT_EQUALS( vtk_mesh_reader2.GetNumNodes(), 12U);
 
#endif //CHASTE_VTK
    }
    void TestMonodomainParallelVtkConversion2D() throw(Exception)
    {
#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        std::string working_directory = "TestHdf5ToVtkConverter_monodomain2D";
        OutputFileHandler handler(working_directory);

        CopyToTestOutputDirectory("heart/test/data/CmguiData/monodomain/2D_0_to_1mm_400_elements.h5",
                                  working_directory);

        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_400_elements");
        DistributedTetrahedralMesh<2,2> mesh(DistributedTetrahedralMeshPartitionType::DUMB);
        mesh.ConstructFromMeshReader(mesh_reader);

        // convert
        HeartConfig::Instance()->SetOutputDirectory(working_directory);
        Hdf5ToVtkConverter<2,2> converter(working_directory, "2D_0_to_1mm_400_elements", &mesh, true);

        //Note that VTK is not thread-safe.  The master process has spawned a child to write the mesh and may still be writing!
        //This barrier just slows things down a bit
        PetscTools::Barrier();

        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        std::stringstream filepath;
        filepath << test_output_directory << working_directory << "/vtk_output/2D_0_to_1mm_400_elements";        
        if (PetscTools::IsSequential())
        {
            filepath <<  ".vtu";
        }
        else
        {
            filepath << "_" << PetscTools::GetMyRank() << ".vtu";
        }
        VtkMeshReader<2,2> vtk_mesh_reader(filepath.str());
        TS_ASSERT_EQUALS( vtk_mesh_reader.GetNumNodes(), mesh.GetNumLocalNodes() + mesh.GetNumHaloNodes() ); //221 in total;
        TS_ASSERT_EQUALS( vtk_mesh_reader.GetNumElements(),  mesh.GetNumLocalElements() ); //400 in total
        
        if (PetscTools::IsSequential())
        {
            std::vector<double> first_node = vtk_mesh_reader.GetNextNode();
            TS_ASSERT_DELTA( first_node[0] , 0.0 , 1e-6 );
            TS_ASSERT_DELTA( first_node[1] , 0.0, 1e-6 );
            TS_ASSERT_DELTA( first_node[2] , 0.0 , 1e-6 );//2d VTK files still carry z-coordinate
    
            std::vector<double> next_node = vtk_mesh_reader.GetNextNode();
            TS_ASSERT_DELTA( next_node[0] , 0.01, 1e-6 );
            TS_ASSERT_DELTA( next_node[1] , 0.0 , 1e-6 );
            TS_ASSERT_DELTA( next_node[2] , 0.0 , 1e-6 );//2d VTK files still carry z-coordinate
        }
        
        //V_m samples 
        std::vector<double> v_at_last;
        vtk_mesh_reader.GetPointData( "V_000020", v_at_last);
        if (PetscTools::IsSequential())
        {
            TS_ASSERT_DELTA( v_at_last[0],   -83.8534, 1e-3 );
            TS_ASSERT_DELTA( v_at_last[110], -83.8534, 1e-3 );
            TS_ASSERT_DELTA( v_at_last[220], -83.8530, 1e-3 );
        }        
        //Show that trying to write pvtu files with original node ordering gives a warning (but writes anyway)
        HeartConfig::Instance()->SetOutputUsingOriginalNodeOrdering(true);
        Hdf5ToVtkConverter<2,2> converter2(working_directory, "2D_0_to_1mm_400_elements", &mesh, true);
        //Note that VTK is not thread-safe.  The master process has spawned a child to write the mesh and may still be writing!
        //This barrier just slows things down a bit
        PetscTools::Barrier();
        VtkMeshReader<2,2> vtk_mesh_reader2(test_output_directory + working_directory +"/vtk_output/2D_0_to_1mm_400_elements.vtu");
        TS_ASSERT_EQUALS( vtk_mesh_reader2.GetNumNodes(), 221U);
  #endif //CHASTE_VTK

    }

    void TestExceptions() throw(Exception)
    {
        std::string directory = "TestHdf5ConverterExceptions";
        HeartConfig::Instance()->SetOutputDirectory(directory);

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements"); //Not used in the test for exceptions until number of nodes is checked
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        CopyToTestOutputDirectory("heart/test/data/CmguiData/monodomain/2D_0_to_1mm_400_elements.h5",
                                  directory);

        TS_ASSERT_THROWS_THIS( CMGUI_3D converter(directory, "2D_0_to_1mm_400_elements", &mesh),
                "Mesh and HDF5 file have a different number of nodes");
    }

};
#endif /*TESTHDF5TOVISUALIZERCONVERTERS_HPP_*/
