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

#include "PostProcessingWriter.hpp"
#include "Hdf5DataReader.hpp"


class TestPropagationPropertiesCalculator : public CxxTest::TestSuite
{
    
public:
    void TestBasicSetup() throw(Exception)
    {
        Hdf5DataReader simulation_data("heart/test/data",
                                       "postprocessingapd", false);
        std::string output_dir = "TestPostProcessingWriter";
        PostProcessingWriter writer(&simulation_data, output_dir);
        
        // These things are available in the XML output requests
        //    <ActionPotentialDurationMap threshold="-30.0" threshold_unit="mV" repolarisation_percentage="90"/>
        //    <UpstrokeTimeMap threshold="-30.0" threshold_unit="mV"/>
        //    <MaxUpstrokeVelocityMap/>
        //    <ConductionVelocityMap origin_node="10"/>
        //    <ConductionVelocityMap origin_node="20"/>
        
        writer.WriteApdMapFile(-30.0, 60.0);
        
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
                                  
        std::string command = "cmp " + test_output_directory + output_dir + "/ApdMap.dat "
                                     + "heart/test/data/good_apd_postprocessing.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);
        
//        writer.WriteUpstrokeTimeMap(-30.0);
//        writer.WriteMaxUpstrokeVelocityMap();
//        writer.WriteConductionVelocityMap(0u);
    }

};


#endif /*TESTPOSTPROCESSINGWRITER_HPP_*/
