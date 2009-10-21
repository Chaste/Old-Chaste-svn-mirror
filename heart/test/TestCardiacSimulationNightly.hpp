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

#ifndef TESTCARDIACSIMULATIONNIGHTLY_HPP_
#define TESTCARDIACSIMULATIONNIGHTLY_HPP_

#include <cxxtest/TestSuite.h>

#include "PetscSetupAndFinalize.hpp"
#include "CardiacSimulation.hpp"

class TestCardiacSimulationNightly : public CxxTest::TestSuite
{
public:
    void TestCardiacSimulationBasicBidomain() throw(Exception)
    {
        // run a bidomain simulation
        CardiacSimulation simulation("heart/test/data/xml/base_bidomain.xml");
        std::string foldername = "BaseBidomain/";
        
        // convert the binary output to a text file using h5dump        
        std::string testoutput_dir = OutputFileHandler::GetChasteTestOutputDirectory();
        std::string command = "h5dump " + testoutput_dir + foldername + "SimulationResults.h5 > " +  testoutput_dir + foldername + "dumped.txt";
        system(command.c_str());
        
        // compare the files, ignoring the first line, which looks like, for example, 'HDF5 "/tmp/chaste/testout/BaseBidomain/SimulationResults.h5" {' 
        command = "diff --ignore-matching-lines=HDF5 heart/test/data/cardiac_simulations/base_bidomain_results.h5 " +  testoutput_dir + foldername + "dumped.txt ";
        int return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value,0); 
    }
    
    void TestCardiacSimulationBasicMonodomain() throw(Exception)
    {
        // run a bidomain simulation
        CardiacSimulation simulation("heart/test/data/xml/base_monodomain.xml");
        std::string foldername = "BaseMonodomain/";
        
        // convert the binary output to a text file using h5dump        
        std::string testoutput_dir = OutputFileHandler::GetChasteTestOutputDirectory();
        std::string command = "h5dump " + testoutput_dir + foldername + "SimulationResults.h5 > " +  testoutput_dir + foldername + "dumped.txt";
        system(command.c_str());
        
        // compare the files, ignoring the first line, which looks like, for example, 'HDF5 "/tmp/chaste/testout/BaseBidomain/SimulationResults.h5" {' 
        command = "diff --ignore-matching-lines=HDF5 heart/test/data/cardiac_simulations/base_monodomain_results.h5 " +  testoutput_dir + foldername + "dumped.txt ";
        int return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value,0); 
    }
    
    void TestCardiacSimulationPostprocessMonodomain() throw(Exception)
    {
        // run a bidomain simulation
        CardiacSimulation simulation("heart/test/data/xml/postprocess_monodomain.xml");
        std::string foldername = "PostprocessMonodomain/";
        
        // convert the binary output to a text file using h5dump        
        std::string testoutput_dir = OutputFileHandler::GetChasteTestOutputDirectory();
        std::string command = "h5dump " + testoutput_dir + foldername + "SimulationResults.h5 > " +  testoutput_dir + foldername + "dumped.txt";
        system(command.c_str());
        
        // compare the files, ignoring the first line, which looks like, for example, 'HDF5 "/tmp/chaste/testout/BaseBidomain/SimulationResults.h5" {' 
        command = "diff --ignore-matching-lines=HDF5 heart/test/data/cardiac_simulations/postprocess_monodomain_results.h5 " +  testoutput_dir + foldername + "dumped.txt ";
        int return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value,0); 
    }
    
    void TestCardiacSimulationSaveBidomain() throw(Exception)
    {
        // run a bidomain simulation
        CardiacSimulation simulation("heart/test/data/xml/save_bidomain.xml");
        std::string foldername = "SaveBidomain";
        
        // convert the binary output to a text file using h5dump        
        std::string testoutput_dir = OutputFileHandler::GetChasteTestOutputDirectory();
        std::string command = "h5dump " + testoutput_dir + foldername + "/SimulationResults.h5 > " +  testoutput_dir + foldername + "/dumped.txt";
        system(command.c_str());
        
        // compare the files, ignoring the first line, which looks like, for example, 'HDF5 "/tmp/chaste/testout/BaseBidomain/SimulationResults.h5" {' 
        command = "diff --ignore-matching-lines=HDF5 heart/test/data/cardiac_simulations/save_bidomain_results.h5 " +  testoutput_dir + foldername + "/dumped.txt ";
        int return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value,0); 

        command = "test -e " +  testoutput_dir + foldername + "_10ms/" + foldername + "_10ms.arch.0"; 
        return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value,0); 
    }
    
    void TestCardiacSimulationSaveMonodomain() throw(Exception)
    {
        // run a bidomain simulation
        CardiacSimulation simulation("heart/test/data/xml/save_monodomain.xml");
        std::string foldername = "SaveMonodomain";
        
        // convert the binary output to a text file using h5dump        
        std::string testoutput_dir = OutputFileHandler::GetChasteTestOutputDirectory();
        std::string command = "h5dump " + testoutput_dir + foldername + "/SimulationResults.h5 > " +  testoutput_dir + foldername + "/dumped.txt";
        system(command.c_str());
        
        // compare the files, ignoring the first line, which looks like, for example, 'HDF5 "/tmp/chaste/testout/BaseBidomain/SimulationResults.h5" {' 
        command = "diff --ignore-matching-lines=HDF5 heart/test/data/cardiac_simulations/save_monodomain_results.h5 " +  testoutput_dir + foldername + "/dumped.txt ";
        int return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value,0); 

        command = "test -e " +  testoutput_dir + foldername + "_10ms/" + foldername + "_10ms.arch.0"; 
        return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value,0); 
    }
};

#endif /*TESTCARDIACSIMULATIONNIGHTLY_HPP_*/
