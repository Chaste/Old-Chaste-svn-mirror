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
    void TestCardiacSimulationNoXmlConstructorNoStim() throw(Exception)
    {
        HeartConfig::Instance()->SetSlabDimensions(0.1, 0.1, 0.1, 0.05);
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-2);
        TS_ASSERT_EQUALS( HeartConfig::Instance()->GetPdeTimeStep(), 0.01);
        HeartConfig::Instance()->SetSimulationDuration(0.01);
        CardiacSimulation simulation;
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "no_stim", false,
                   "ChasteResults", "SimulationResults", true));
        HeartConfig::Instance()->Reset();
    }

    void TestMono1dNoStim() throw(Exception)
    {
        HeartConfig::Instance()->SetFibreLength(0.1, 0.05);
        HeartConfig::Instance()->SetSimulationDuration(0.01);
        HeartConfig::Instance()->SetSpaceDimension(1);
        HeartConfig::Instance()->SetDomain(cp::domain_type::Mono);
        CardiacSimulation simulation;
        HeartConfig::Instance()->Reset();
    }    
    void TestMono2dNoStim() throw(Exception)
    {
        HeartConfig::Instance()->SetSheetDimensions(0.1, 0.1, 0.05);
        HeartConfig::Instance()->SetSimulationDuration(0.01);
        HeartConfig::Instance()->SetSpaceDimension(2);
        HeartConfig::Instance()->SetDomain(cp::domain_type::Mono);
        CardiacSimulation simulation;
        HeartConfig::Instance()->Reset();
    }    
    void TestMono3dNoStim() throw(Exception)
    {
        HeartConfig::Instance()->SetSlabDimensions(0.1, 0.1, 0.1, 0.05);
        HeartConfig::Instance()->SetSimulationDuration(0.01);
        HeartConfig::Instance()->SetSpaceDimension(3);
        HeartConfig::Instance()->SetDomain(cp::domain_type::Mono);
        HeartConfig::Instance()->SetOutputDirectory("Mono3dNoStim");
        HeartConfig::Instance()->SetSaveSimulation(true);
        CardiacSimulation simulation;
        HeartConfig::Instance()->Reset();
    }    

    void TestBi1dNoStim() throw(Exception)
    {
        HeartConfig::Instance()->SetFibreLength(0.1, 0.05);
        HeartConfig::Instance()->SetSimulationDuration(0.01);
        HeartConfig::Instance()->SetSpaceDimension(1);
        HeartConfig::Instance()->SetDomain(cp::domain_type::Bi);
        CardiacSimulation simulation;
        HeartConfig::Instance()->Reset();
    }    
    void TestBi2dNoStim() throw(Exception)
    {
        HeartConfig::Instance()->SetSheetDimensions(0.1, 0.1, 0.05);
        HeartConfig::Instance()->SetSimulationDuration(0.01);
        HeartConfig::Instance()->SetSpaceDimension(2);
        HeartConfig::Instance()->SetDomain(cp::domain_type::Bi);
        CardiacSimulation simulation;
        HeartConfig::Instance()->Reset();
    }    
    void TestBi3dNoStim() throw(Exception)
    {
        HeartConfig::Instance()->SetSlabDimensions(0.1, 0.1, 0.1, 0.05);
        HeartConfig::Instance()->SetSimulationDuration(0.01);
        HeartConfig::Instance()->SetSpaceDimension(3);
        HeartConfig::Instance()->SetDomain(cp::domain_type::Bi);
        CardiacSimulation simulation;
        HeartConfig::Instance()->Reset();
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
   void TestCardiacSimulationPatchwork() throw(Exception)
    {
        // run a bidomain simulation
        CardiacSimulation simulation("heart/test/data/xml/base_monodomain_patchwork.xml");
        std::string foldername = "Patchwork";
        
       // compare the files, using the CompareFilesViaHdf5DataReader() method  
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "patchwork_results", false,
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
    
    void TestCardiacSimulationSaveBidomain() throw(Exception)
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
    
    void TestCardiacSimulationSaveMonodomain() throw(Exception)
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
    
    void TestExceptions() throw(Exception)
    {
        HeartConfig::Instance()->SetSpaceDimension(8); //Electro-physio-super-string
        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation, "Bidomain space dimension not supported: should be 1, 2 or 3");
        HeartConfig::Instance()->SetDomain(cp::domain_type::Mono);
        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation, "Monodomain space dimension not supported: should be 1, 2 or 3");
        
        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation("no file"),"Missing file parsing configuration file: no file");
    }
};

#endif /*TESTCARDIACSIMULATION_HPP_*/

