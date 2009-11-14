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

#ifndef TESTCARDIACSIMULATION_HPP_
#define TESTCARDIACSIMULATION_HPP_

#include <cxxtest/TestSuite.h>

#include "CardiacSimulation.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "CompareHdf5ResultsFiles.hpp"

class TestCardiacSimulation : public CxxTest::TestSuite
{
public:

    void TestMono1dSmall() throw(Exception)
    {
        CardiacSimulation simulation("heart/test/data/xml/monodomain1d_small.xml");
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "mono_1d_small", false,
                   "SaveMono1D", "SimulationResults", true));
        CardiacSimulation simulation2("heart/test/data/xml/monodomain1d_resume.xml");
    }
    void TestMono2dSmall() throw(Exception)
    {
        CardiacSimulation simulation("heart/test/data/xml/monodomain2d_small.xml");
        CardiacSimulation simulation2("heart/test/data/xml/monodomain2d_resume.xml");
    }
    void TestMono3dSmall() throw(Exception)
    {
        CardiacSimulation simulation("heart/test/data/xml/monodomain3d_small.xml");
        CardiacSimulation simulation2("heart/test/data/xml/monodomain3d_resume.xml");
    }
    void TestBi1dSmall() throw(Exception)
    {
        CardiacSimulation simulation("heart/test/data/xml/bidomain1d_small.xml");
        CardiacSimulation simulation2("heart/test/data/xml/bidomain1d_resume.xml");
    }
    void TestBi2dSmall() throw(Exception)
    {
        CardiacSimulation simulation("heart/test/data/xml/bidomain2d_small.xml");
        CardiacSimulation simulation2("heart/test/data/xml/bidomain2d_resume.xml");
    }
    void TestBi3dSmall() throw(Exception)
    {
        CardiacSimulation simulation("heart/test/data/xml/bidomain3d_small.xml");
        CardiacSimulation simulation2("heart/test/data/xml/bidomain3d_resume.xml");
    }

    void TestCardiacSimulationBasicBidomainShort() throw(Exception)
    {
        // run a bidomain simulation
        CardiacSimulation simulation("heart/test/data/xml/base_bidomain_short.xml");

        // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "base_bidomain_short_results", false,
                   "BaseBidomainShort", "SimulationResults", true));
    }
    
    void TestCardiacSimulationBasicMonodomainShort() throw(Exception)
    {
        // run a bidomain simulation
        CardiacSimulation simulation("heart/test/data/xml/base_monodomain_short.xml");
        std::string foldername = "BaseMonodomainShort";
        
       // compare the files, using the CompareFilesViaHdf5DataReader() method  
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "base_monodomain_short_results", false,
                   foldername, "SimulationResults", true));
    }
    
    void TestCardiacSimulationPostprocessMonodomain() throw(Exception)
    {
        // run a bidomain simulation
        CardiacSimulation simulation("heart/test/data/xml/postprocess_monodomain_short.xml");
        std::string foldername = "PostprocessMonodomainShort";
        
        // compare the files, using the CompareFilesViaHdf5DataReader() method  
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "postprocess_monodomain_short_results", false,
                   foldername, "SimulationResults", true));
    }
    
    void TestCardiacSimulationArchiveBidomain() throw(Exception)
    {
        // run a bidomain simulation
        CardiacSimulation simulation("heart/test/data/xml/save_bidomain_short.xml");
        std::string foldername = "SaveBidomainShort";
        
        // compare the files, using the CompareFilesViaHdf5DataReader() method  
        TS_ASSERT(CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "save_bidomain_short_results", false,
                                                foldername, "SimulationResults", true));

        std::string command = "test -e " +  OutputFileHandler::GetChasteTestOutputDirectory() + foldername + "_0.2ms/" + foldername + "_0.2ms.arch.0"; 
        int return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value,0); 
    }
    
    void TestCardiacSimulationArchiveMonodomain() throw(Exception)
    {
        // run a bidomain simulation
        CardiacSimulation simulation("heart/test/data/xml/save_monodomain_short.xml");
        std::string foldername = "SaveMonodomainShort";
        
        // compare the files, using the CompareFilesViaHdf5DataReader() method  
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "save_monodomain_short_results", false,
                   foldername, "SimulationResults", true));

        std::string command = "test -e " +  OutputFileHandler::GetChasteTestOutputDirectory() + foldername + "_0.2ms/" + foldername + "_0.2ms.arch.0"; 
        int return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value,0); 
    }
    
    // requires TestCardiacSimulationSaveMonodomain() to have been run
    void TestCardiacSimulationResumeMonodomain() throw(Exception)
    {
        // run a bidomain simulation
        HeartConfig::Instance()->SetSpaceDimension(1);
        CardiacSimulation simulation("heart/test/data/xml/resume_monodomain_short.xml");
        std::string foldername = "SaveMonodomainShort";
        
        // compare the files, using the CompareFilesViaHdf5DataReader() method  
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "resume_monodomain_short_results", false,
                   foldername, "SimulationResults", true));
    }
    
    // requires TestCardiacSimulationSaveBidomain() to have been run
    void TestCardiacSimulationResumeBidomain() throw(Exception)
    {
        // run a bidomain simulation
        HeartConfig::Instance()->SetSpaceDimension(1);
        CardiacSimulation simulation("heart/test/data/xml/resume_bidomain_short.xml");
        std::string foldername = "SaveBidomainShort";
        
        // compare the files, using the CompareFilesViaHdf5DataReader() method  
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "resume_bidomain_short_results", false,
                   foldername, "SimulationResults", true));
    }
    
   void TestCardiacSimulationPatchwork() throw(Exception)
    {
        CardiacSimulation simulation("heart/test/data/xml/base_monodomain_patchwork.xml");
        std::string foldername = "Patchwork";
        
       // compare the files, using the CompareFilesViaHdf5DataReader() method  
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "patchwork_results", false,
                   foldername, "SimulationResults", true));
    }
    void TestCardiacSimulationKirsten() throw(Exception)
    {
        CardiacSimulation simulation("heart/test/data/xml/base_monodomain_kirsten.xml");
    }

    void TestExceptions() throw(Exception)
    {
        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation("heart/test/data/xml/monodomain8d_small.xml"), "Monodomain space dimension not supported: should be 1, 2 or 3");
        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation("heart/test/data/xml/bidomain8d_small.xml"), "Bidomain space dimension not supported: should be 1, 2 or 3");

        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation("heart/test/data/xml/base_monodomain_frankenstein.xml"),
          "XML parsing error in configuration file: heart/test/data/xml/base_monodomain_frankenstein.xml");
        
        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation("no file"),"Missing file parsing configuration file: no file");
        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation(""),"No XML file name given");
    }
};

#endif /*TESTCARDIACSIMULATION_HPP_*/

