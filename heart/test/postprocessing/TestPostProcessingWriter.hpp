/*

Copyright (c) 2005-2013, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef TESTPOSTPROCESSINGWRITER_HPP_
#define TESTPOSTPROCESSINGWRITER_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>

#include "PostProcessingWriter.hpp"
#include "OutputFileHandler.hpp"
#include "FileFinder.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "DistanceMapCalculator.hpp"
#include "HeartConfig.hpp"
#include "NumericFileComparison.hpp"
#include "BidomainProblem.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "PetscSetupAndFinalize.hpp"

//#include "VtkMeshReader.hpp" //Needed for commented out test, see #1660


class TestPostProcessingWriter : public CxxTest::TestSuite
{
    /**
     * Set the output directory in HeartConfig, and return where our files of interest are.
     * @param rOutputDirName  the leaf folder name
     * @return  FileFinder for the output folder within this
     */
    FileFinder GetOutputPath(const std::string& rOutputDirName)
    {
        HeartConfig::Instance()->SetOutputDirectory(rOutputDirName);
        OutputFileHandler handler(rOutputDirName);
        return handler.FindFile("output");
    }

public:
    void TestWriterMethods() throw(Exception)
    {
        FileFinder output_dir = GetOutputPath("TestPostProcessingWriter_WriterMethods");

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        PostProcessingWriter<1,1> writer(mesh, "heart/test/data", "postprocessingapd", false);

        writer.WriteApdMapFile(60.0, -30.0);

        std::string file1 = FileFinder("Apd_60_-30_Map.dat", output_dir).GetAbsolutePath();
        std::string file2 = "heart/test/data/PostProcessorWriter/good_apd_postprocessing.dat";
        NumericFileComparison comp(file1, file2);
        TS_ASSERT(comp.CompareFiles(1e-12));

        writer.WriteUpstrokeTimeMap(-30.0);

        file1 = FileFinder("UpstrokeTimeMap_-30.dat", output_dir).GetAbsolutePath();
        file2 = "heart/test/data/PostProcessorWriter/good_upstroke_time_postprocessing.dat";
        NumericFileComparison comp2(file1, file2);
        TS_ASSERT(comp2.CompareFiles(1e-12));

        writer.WriteMaxUpstrokeVelocityMap(-30.0);

        file1 = FileFinder("MaxUpstrokeVelocityMap_-30.dat", output_dir).GetAbsolutePath();
        file2 = "heart/test/data/PostProcessorWriter/good_upstroke_velocity_postprocessing.dat";
        NumericFileComparison comp3(file1, file2);
        TS_ASSERT(comp3.CompareFiles(1e-12));

        DistanceMapCalculator<1,1> dist_calculator(mesh);

        std::vector<unsigned> origin_node;
        origin_node.push_back(0);
        std::vector<double> distance_map_from_0;

        dist_calculator.ComputeDistanceMap(origin_node, distance_map_from_0);

        writer.WriteConductionVelocityMap(0u, distance_map_from_0);

        file1 = FileFinder("ConductionVelocityFromNode0.dat", output_dir).GetAbsolutePath();
        file2 = "heart/test/data/PostProcessorWriter/conduction_velocity_10_nodes_from_node_0.dat";
        NumericFileComparison comp4(file1, file2);
        TS_ASSERT(comp4.CompareFiles(1e-12));
    }

    void TestApdWritingWithNoApdsPresent() throw(Exception)
    {
        FileFinder output_dir = GetOutputPath("TestPostProcessingWriter_ApdWritingWithNoApdsPresent");

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_100_elements");
        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        PostProcessingWriter<1,1> writer(mesh, "heart/test/data/Monodomain1d", "MonodomainLR91_1d", false);

        writer.WriteApdMapFile(90.0, -30.0);

        std::string file1 = FileFinder("Apd_90_-30_Map.dat", output_dir).GetAbsolutePath();
        std::string file2 = "heart/test/data/PostProcessorWriter/101_zeroes.dat";
        NumericFileComparison comp(file1, file2);
        TS_ASSERT(comp.CompareFiles(1e-12));
    }

    void TestPostProcessWriting() throw (Exception)
    {
        FileFinder output_dir = GetOutputPath("TestPostProcessingWriter_PostProcessWriting");

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_10_100_elements");
        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<std::pair<double,double> > apd_maps;
        apd_maps.push_back(std::pair<double, double>(80,-30));//repolarisation percentage first, as per schema
        apd_maps.push_back(std::pair<double, double>(90,-20));//repolarisation percentage first, as per schema
        HeartConfig::Instance()->SetApdMaps(apd_maps);

        std::vector<double> upstroke_time_map;
        upstroke_time_map.push_back(-70.0);
        upstroke_time_map.push_back( 20.0);
        HeartConfig::Instance()->SetUpstrokeTimeMaps(upstroke_time_map);

        std::vector<double> upstroke_velocity_map;
        upstroke_velocity_map.push_back(-50.0);
        upstroke_velocity_map.push_back(50.0);
        HeartConfig::Instance()->SetMaxUpstrokeVelocityMaps(upstroke_velocity_map);

        std::vector<unsigned> conduction_velocity_map;
        conduction_velocity_map.push_back(0u);
        HeartConfig::Instance()->SetConductionVelocityMaps(conduction_velocity_map);

        //test the mtehod that extrapolates nodal traces
        std::vector<unsigned> nodes_to_extrapolate;//1D test, does not cover node permutation case
        nodes_to_extrapolate.push_back(1u);
        nodes_to_extrapolate.push_back(99u);
        HeartConfig::Instance()->SetRequestedNodalTimeTraces(nodes_to_extrapolate);

        std::vector<ChastePoint<1> > pseudo_ecg_electrodes;
        pseudo_ecg_electrodes.push_back(ChastePoint<1>(11.0));
        pseudo_ecg_electrodes.push_back(ChastePoint<1>(-1.0));
        HeartConfig::Instance()->SetPseudoEcgElectrodePositions(pseudo_ecg_electrodes);

        PostProcessingWriter<1,1> writer(mesh, "heart/test/data/Monodomain1d", "MonodomainLR91_1d", false);

        writer.WritePostProcessingFiles();

        writer.WriteAboveThresholdDepolarisationFile(-40.0);

        std::string file1 = FileFinder("Apd_80_-30_Map.dat", output_dir).GetAbsolutePath();
        std::string file2 = "heart/test/data/PostProcessorWriter/101_zeroes.dat";
        NumericFileComparison comp1(file1, file2);
        TS_ASSERT(comp1.CompareFiles(1e-12));

        file1 = FileFinder("Apd_90_-20_Map.dat", output_dir).GetAbsolutePath();
        file2 = "heart/test/data/PostProcessorWriter/101_zeroes.dat";
        NumericFileComparison comp2(file1, file2);
        TS_ASSERT(comp2.CompareFiles(1e-12));

        file1 = FileFinder("UpstrokeTimeMap_-70.dat", output_dir).GetAbsolutePath();
        file2 = "heart/test/data/PostProcessorWriter/UpstrokeTimeMap_-70.dat";
        NumericFileComparison comp3(file1, file2);
        TS_ASSERT(comp3.CompareFiles(1e-12));

        file1 = FileFinder("UpstrokeTimeMap_20.dat", output_dir).GetAbsolutePath();
        file2 = "heart/test/data/PostProcessorWriter/UpstrokeTimeMap_20.dat";
        NumericFileComparison comp4(file1, file2);
        TS_ASSERT(comp4.CompareFiles(1e-12));

        file1 = FileFinder("MaxUpstrokeVelocityMap_-50.dat", output_dir).GetAbsolutePath();
        file2 = "heart/test/data/PostProcessorWriter/MaxUpstrokeVelocityMap_-50.dat";
        NumericFileComparison comp5(file1, file2);
        TS_ASSERT(comp5.CompareFiles(1e-12));

        file1 = FileFinder("MaxUpstrokeVelocityMap_50.dat", output_dir).GetAbsolutePath();
        file2 = "heart/test/data/PostProcessorWriter/MaxUpstrokeVelocityMap_50.dat";
        NumericFileComparison comp6(file1, file2);
        TS_ASSERT(comp6.CompareFiles(1e-12));

        file1 = FileFinder("ConductionVelocityFromNode0.dat", output_dir).GetAbsolutePath();
        file2 = "heart/test/data/PostProcessorWriter/conduction_velocity_100_nodes_from_node_0.dat";
        NumericFileComparison comp7(file1, file2);
        TS_ASSERT(comp7.CompareFiles(1e-12));

        file1 = FileFinder("PseudoEcgFromElectrodeAt_11_0_0.dat", output_dir).GetAbsolutePath();
        file2 = "heart/test/data/PostProcessorWriter/PseudoEcgFromElectrodeAt_11_0_0.dat";
        NumericFileComparison comp8(file1, file2);
        TS_ASSERT(comp8.CompareFiles(1e-12));

        file1 = FileFinder("PseudoEcgFromElectrodeAt_-1_0_0.dat", output_dir).GetAbsolutePath();
        file2 = "heart/test/data/PostProcessorWriter/PseudoEcgFromElectrodeAt_-1_0_0.dat";
        NumericFileComparison comp9(file1, file2);
        TS_ASSERT(comp9.CompareFiles(1e-12));

        file1 = FileFinder("AboveThresholdDepolarisations-40.dat", output_dir).GetAbsolutePath();
        file2 = "heart/test/data/PostProcessorWriter/AboveThresholdDepolarisations-40.dat";
        NumericFileComparison comp10(file1, file2);
        TS_ASSERT(comp10.CompareFiles(1e-12));

        file1 = FileFinder("NodalTraces_V.dat", output_dir).GetAbsolutePath();
        file2 = "heart/test/data/PostProcessorWriter/NodalTrace_V_Valid.dat";
        NumericFileComparison comp_nodes(file1, file2);
        TS_ASSERT(comp_nodes.CompareFiles(1e-12));
    }

    void TestExtractNodeTracesWithNodePermutation() throw (Exception)
    {
        HeartConfig::Instance()->Reset();
        FileFinder output_dir = GetOutputPath("TestPostProcessingWriter_ExtractNodeTracesWithNodePermutation");
        HeartConfig::Instance()->SetOutputFilenamePrefix("NodalTracesTest");

        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_400_elements");
        DistributedTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // The point of this test is to check the method handles permutation properly...
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetOutputUsingOriginalNodeOrdering(), false);
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005, 0.0005));
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);
        HeartConfig::Instance()->SetSimulationDuration(2); //ms

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory;
        BidomainProblem<2> problem( &cell_factory );
        problem.SetMesh(&mesh);

        std::vector<unsigned> nodes;
        nodes.push_back(0u);
        nodes.push_back(13u);
        HeartConfig::Instance()->SetRequestedNodalTimeTraces(nodes);

        problem.Initialise();
        problem.Solve();

        // compare for two nodes (0 and 13) that I visually checked after running a no-permutation simulation (with Tetrahedral mesh)
        // and comparing the hdf file (with hdfview) with the postprocessing output
        // in order to consider the file "valid".

        //check the file with V (assuming default variable name was V)
        std::string file1 = FileFinder("NodalTraces_V.dat", output_dir).GetAbsolutePath();
        std::string file2 = "heart/test/data/PostProcessorWriter/NodalTrace_V_WithPermutationValid.dat";

        NumericFileComparison comp_V(file1, file2);
        TS_ASSERT(comp_V.CompareFiles(1e-3));

        //check file with Phi_e (assuming default variable name was Phi_e)
        file1 = FileFinder("NodalTraces_Phi_e.dat", output_dir).GetAbsolutePath();
        file2 = "heart/test/data/PostProcessorWriter/NodalTrace_Phi_e_WithPermutationValid.dat";

        NumericFileComparison comp_phie(file1, file2);
        TS_ASSERT(comp_phie.CompareFiles(1e-3));
    }

    void TestWritingEads() throw (Exception)
    {
        HeartConfig::Instance()->Reset();
        FileFinder output_dir = GetOutputPath("TestPostProcessingWriter_WritingEads");

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_10_100_elements");
        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        PostProcessingWriter<1,1> writer(mesh, "heart/test/data", "Ead", false);

        writer.WriteAboveThresholdDepolarisationFile(-30.0);

        std::string file1 = FileFinder("AboveThresholdDepolarisations-30.dat", output_dir).GetAbsolutePath();
        std::string file2 = "heart/test/data/PostProcessorWriter/AboveThresholdDepolarisations-30.dat";
        NumericFileComparison comp(file1, file2);
        TS_ASSERT(comp.CompareFiles(1e-12));
    }

    void TestSwitchingOutputFormat() throw (Exception)
    {
        HeartConfig::Instance()->Reset();
        FileFinder output_dir = GetOutputPath("TestPostProcessingWriter_SwitchingOutputFormat");
        HeartConfig::Instance()->SetVisualizeWithCmgui();

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        PostProcessingWriter<1,1> writer(mesh, "heart/test/data", "postprocessingapd", false);

        writer.WriteApdMapFile(60.0, -30.0);

        std::string file1 = FileFinder("Apd_60_-30_Map.dat", output_dir).GetAbsolutePath();
        std::string file2 = "heart/test/data/PostProcessorWriter/good_apd_postprocessing.dat";
        NumericFileComparison comp(file1, file2);
        TS_ASSERT(comp.CompareFiles(1e-12));
    }

// Test fails as VTK post processing output not yet implemented - See #1660
//    void xxxTestVtkOutput() throw (Exception)
//    {
//#ifdef CHASTE_VTK
//        HeartConfig::Instance()->Reset();
//
//        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_400_elements");
//        DistributedTetrahedralMesh<2,2> mesh;
//        mesh.ConstructFromMeshReader(mesh_reader);
//
//        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005, 0.0005));
//        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
//        HeartConfig::Instance()->SetCapacitance(1.0);
//        HeartConfig::Instance()->SetSimulationDuration(2); //ms
//        HeartConfig::Instance()->SetVisualizeWithVtk();
//        HeartConfig::Instance()->SetOutputDirectory("TestPostProcessingWriter_VtkOutput");
//
//        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory;
//        BidomainProblem<2> problem( &cell_factory );
//        problem.SetMesh(&mesh);
//
//        problem.Initialise();
//        problem.Solve();
//
//        OutputFileHandler handler("TestPostProcessingWriter_VtkOutput", false);
//        FileFinder output_dir = handler.FindFile("vtk_output");
//
//        PostProcessingWriter<2,2> writer(mesh, handler.GetOutputDirectoryFullPath(), HeartConfig::Instance()->GetOutputFilenamePrefix(), false);
//
//        writer.WriteApdMapFile(60.0, -30.0);
//
//        //Read VTK file & check data.
//        VtkMeshReader<2,2> vtk_reader(output_dir.GetAbsolutePath() + "/SimulationResults.vtu");
//
//        std::vector<double> apd_data;
//
//        //This will currently fail as the APD map data is not yet written to the VTK file.
//        vtk_reader.GetPointData("ApdMap", apd_data);
//#endif
//    }
};


#endif /*TESTPOSTPROCESSINGWRITER_HPP_*/
