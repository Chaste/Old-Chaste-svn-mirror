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

#ifndef TESTPOSTPROCESSINGWRITER_HPP_
#define TESTPOSTPROCESSINGWRITER_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>

#include "NonCachedTetrahedralMesh.hpp" //must be first, it gets UblasIncludes from the mesh classes (ChastePoint.hpp)
#include "PostProcessingWriter.hpp"
#include "Hdf5DataReader.hpp"
#include "OutputFileHandler.hpp"
#include "TrianglesMeshReader.hpp"
#include "DistanceMapCalculator.hpp"
#include "HeartConfig.hpp"
#include "TrianglesMeshReader.hpp"
#include "ParallelTetrahedralMesh.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestPostProcessingWriter : public CxxTest::TestSuite
{
    
public:
// These things are available in the XML output requests
//    <ActionPotentialDurationMap threshold="-30.0" threshold_unit="mV" repolarisation_percentage="90"/>
//    <UpstrokeTimeMap threshold="-30.0" threshold_unit="mV"/>
//    <MaxUpstrokeVelocityMap/>
//    <ConductionVelocityMap origin_node="10"/>
//    <ConductionVelocityMap origin_node="20"/>
        
    void TestWriterMethods() throw(Exception)
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements"); 
        ParallelTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        std::string output_dir = "ChasteResults/output"; // default given by HeartConfig
        PostProcessingWriter<1,1> writer(mesh, "heart/test/data", "postprocessingapd", false);
        
        writer.WriteApdMapFile(60.0, -30.0);
                                   
        std::string command;
        command = "cmp " + OutputFileHandler::GetChasteTestOutputDirectory() 
                                    + output_dir + "/Apd_60_-30_Map.dat "
                                    + "heart/test/data/PostProcessorWriter/good_apd_postprocessing.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);
        
        writer.WriteUpstrokeTimeMap(-30.0);
        command = "cmp " + OutputFileHandler::GetChasteTestOutputDirectory() 
                                    + output_dir + "/UpstrokeTimeMap_-30.dat "
                                    + "heart/test/data/PostProcessorWriter/good_upstroke_time_postprocessing.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);
        
        writer.WriteMaxUpstrokeVelocityMap(-30.0);
        command = "cmp " + OutputFileHandler::GetChasteTestOutputDirectory() 
                                    + output_dir + "/MaxUpstrokeVelocityMap_-30.dat "
                                    + "heart/test/data/PostProcessorWriter/good_upstroke_velocity_postprocessing.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);

        DistanceMapCalculator<1,1> dist_calculator(mesh);
        
        std::vector<unsigned> origin_node;
        origin_node.push_back(0);
        std::vector<double> distance_map_from_0;
        
        dist_calculator.ComputeDistanceMap(origin_node, distance_map_from_0);

        writer.WriteConductionVelocityMap(0u, distance_map_from_0);
        command = "cmp " + OutputFileHandler::GetChasteTestOutputDirectory() 
                                    + output_dir + "/ConductionVelocityFromNode0.dat "
                                    + "heart/test/data/PostProcessorWriter/conduction_velocity_10_nodes_from_node_0.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);
    }
    
    void TestApdWritingWithNoApdsPresent() throw(Exception)
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_100_elements"); 
        ParallelTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        std::string output_dir = "ChasteResults/output"; // default given by HeartConfig
        PostProcessingWriter<1,1> writer(mesh, "heart/test/data/Monodomain1d", "MonodomainLR91_1d", false);
        
        writer.WriteApdMapFile(90.0, -30.0);
        
        std::string command = "cmp " + OutputFileHandler::GetChasteTestOutputDirectory() + output_dir + "/Apd_90_-30_Map.dat " 
                                   + "heart/test/data/101_zeroes.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);
        
    }
    
    void TestPostProcessWriting() throw (Exception)
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_10_100_elements"); 
        ParallelTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        std::vector<std::pair<double,double> > apd_maps;
        apd_maps.push_back(std::pair<double, double>(80,-30));//reploarisation percentage first, as per schema
        apd_maps.push_back(std::pair<double, double>(90,-20));//reploarisation percentage first, as per schema
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
                                                                            
        PostProcessingWriter<1,1> writer(mesh, "heart/test/data/Monodomain1d", "MonodomainLR91_1d", false);  

        writer.WritePostProcessingFiles();
        
        std::string output_dir = "ChasteResults/output"; // default given by HeartConfig
        std::string command = "cmp " + OutputFileHandler::GetChasteTestOutputDirectory() + output_dir + "/Apd_80_-30_Map.dat " 
                                   + "heart/test/data/101_zeroes.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);                         

        command = "cmp " + OutputFileHandler::GetChasteTestOutputDirectory() + output_dir + "/Apd_90_-20_Map.dat " 
                  + "heart/test/data/101_zeroes.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);      

        command = "cmp " + OutputFileHandler::GetChasteTestOutputDirectory() + output_dir + "/UpstrokeTimeMap_-70.dat " 
                  + "heart/test/data/PostProcessorWriter/UpstrokeTimeMap_-70.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);

        command = "cmp " + OutputFileHandler::GetChasteTestOutputDirectory() + output_dir + "/UpstrokeTimeMap_20.dat " 
                  + "heart/test/data/PostProcessorWriter/UpstrokeTimeMap_20.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);

        command = "cmp " + OutputFileHandler::GetChasteTestOutputDirectory() + output_dir + "/MaxUpstrokeVelocityMap_-50.dat " 
                  + "heart/test/data/PostProcessorWriter/MaxUpstrokeVelocityMap_-50.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);

        command = "cmp " + OutputFileHandler::GetChasteTestOutputDirectory() + output_dir + "/MaxUpstrokeVelocityMap_50.dat " 
                  + "heart/test/data/PostProcessorWriter/MaxUpstrokeVelocityMap_50.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);

        command = "cmp " + OutputFileHandler::GetChasteTestOutputDirectory() + output_dir + "/ConductionVelocityFromNode0.dat " 
                  + "heart/test/data/PostProcessorWriter/conduction_velocity_100_nodes_from_node_0.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);
        
        
    }
};


#endif /*TESTPOSTPROCESSINGWRITER_HPP_*/
