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

#ifndef TESTBIDOMAINARCHIVEKSP_HPP_
#define TESTBIDOMAINARCHIVEKSP_HPP_


#include <cxxtest/TestSuite.h>


/// \todo #98: test unarchiving
//#include <boost/archive/text_iarchive.hpp>

#include "LinearSystem.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "BidomainProblem.hpp"
#include "PetscSetupAndFinalize.hpp"


class TestBidomainArchiveLinearSystem : public CxxTest::TestSuite
{
public:

    /**
     *  \todo #98: test unarchiving
     */
    void TestDumpBidomainPE() throw (Exception)
    {
        // after 1.4 ms of simulation the wavefront is half way through the cube
        HeartConfig::Instance()->SetSimulationDuration(1.4);  //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/3D_0_to_1mm_6000_elements");
        HeartConfig::Instance()->SetOutputDirectory("DumpBidomain3DPE");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomainPE3d");
        HeartConfig::Instance()->SetVisualizeWithMeshalyzer(true);

        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 3> bidomain_cell_factory(-600.0*1000);

        BidomainProblem<3> bidomain_problem( &bidomain_cell_factory );

        bidomain_problem.SetArchiveLinearSystemObject();   
        bidomain_problem.Initialise();
        bidomain_problem.Solve();
        
//        OutputFileHandler handler(HeartConfig::Instance()->GetOutputDirectory(), false);
//        handler.SetArchiveDirectory();
//        
//        std::string archive_filename;
//        archive_filename = handler.GetOutputDirectoryFullPath() + HeartConfig::Instance()->GetOutputFilenamePrefix() + "_ls.arch";       
//
//        std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
//        boost::archive::text_iarchive input_arch(ifs); 
//        
//        //LinearSystem linear_system(3);
//        LinearSystem* p_linear_system;
//        input_arch >> p_linear_system;
//        
//        delete p_linear_system;
    }
    
    
};

#endif /*TESTBIDOMAINARCHIVEKSP_HPP_*/
