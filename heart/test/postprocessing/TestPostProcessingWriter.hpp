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

class TestPostProcessingWriter : public CxxTest::TestSuite
{
    
public:
// These things are available in the XML output requests
//    <ActionPotentialDurationMap threshold="-30.0" threshold_unit="mV" repolarisation_percentage="90"/>
//    <UpstrokeTimeMap threshold="-30.0" threshold_unit="mV"/>
//    <MaxUpstrokeVelocityMap/>
//    <ConductionVelocityMap origin_node="10"/>
//    <ConductionVelocityMap origin_node="20"/>
        
    void TestWriting() throw(Exception)
    {
        Hdf5DataReader simulation_data("heart/test/data",
                                       "postprocessingapd", false);
                                       
        std::string output_dir = "ChasteResults/output"; // default given by HeartConfig
        PostProcessingWriter writer(&simulation_data);
        
        writer.WriteApdMapFile(-30.0, 60.0);
                                   
        std::string command;
        command = "cmp " + OutputFileHandler::GetChasteTestOutputDirectory() 
                                    + output_dir + "/Apd60Map.dat "
                                    + "heart/test/data/good_apd_postprocessing.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);
        
        writer.WriteUpstrokeTimeMap(-30.0);
        command = "cmp " + OutputFileHandler::GetChasteTestOutputDirectory() 
                                    + output_dir + "/UpstrokeTimeMap.dat "
                                    + "heart/test/data/good_upstroke_time_postprocessing.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);
        
        writer.WriteMaxUpstrokeVelocityMap(-30.0);
        command = "cmp " + OutputFileHandler::GetChasteTestOutputDirectory() 
                                    + output_dir + "/MaxUpstrokeVelocityMap.dat "
                                    + "heart/test/data/good_upstroke_velocity_postprocessing.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        NonCachedTetrahedralMesh<1,1> cable_mesh;
        cable_mesh.ConstructFromMeshReader(mesh_reader);

        DistanceMapCalculator<1> dist_calculator(cable_mesh);
        
        std::vector<unsigned> origin_node;
        origin_node.push_back(0);
        std::vector<double> distance_map_from_0;
        
        dist_calculator.ComputeDistanceMap(origin_node, distance_map_from_0);

        writer.WriteConductionVelocityMap(0u, distance_map_from_0);
        command = "cmp " + OutputFileHandler::GetChasteTestOutputDirectory() 
                                    + output_dir + "/ConductionVelocityFromNode0.dat "
                                    + "heart/test/data/good_conduction_velocity_postprocessing.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);
    }
    
    void TestApdWritingWithNoApdsPresent() throw(Exception)
    {
        Hdf5DataReader simulation_data("heart/test/data/Monodomain1d",
                                       "MonodomainLR91_1d", false);
                                       
        std::string output_dir = "ChasteResults/output"; // default given by HeartConfig
        PostProcessingWriter writer(&simulation_data);  
        
        writer.WriteApdMapFile(-30.0, 90.0);
        
        std::string command = "cmp " + OutputFileHandler::GetChasteTestOutputDirectory() + output_dir + "/Apd90Map.dat " 
                                   + "heart/test/data/101_zeroes.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);
    }
    
    void xTestPostProcessWriting() throw (Exception)
    {
        //HeartConfig::Instance()->Set...
        
        //Constructor writes info based on what's in the config
        Hdf5DataReader simulation_data("heart/test/data/Monodomain1d",
                                       "MonodomainLR91_1d", false);
        std::string output_dir = "ChasteResults/output"; // default given by HeartConfig
        PostProcessingWriter writer(&simulation_data);  
        
            std::string command = "cmp " + OutputFileHandler::GetChasteTestOutputDirectory() + output_dir + "/Apd90Map.dat " 
                                   + "heart/test/data/101_zeroes.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);                         
                                       
    }
    
    

};


#endif /*TESTPOSTPROCESSINGWRITER_HPP_*/
