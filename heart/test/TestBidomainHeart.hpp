/*

Copyright (C) University of Oxford, 2008

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


#ifndef _TESTBIDOMAINHEART_HPP_
#define _TESTBIDOMAINHEART_HPP_

#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include <petscvec.h>
#include <vector>

#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "TrianglesMeshWriter.hpp"
#include "PetscTools.hpp"



class PointStimulusHeartCellFactory : public AbstractCardiacCellFactory<3>
{
private:
    SimpleStimulus *mpStimulus;
public:
    PointStimulusHeartCellFactory() : AbstractCardiacCellFactory<3>()
    {
        mpStimulus = new SimpleStimulus(-1000.0*500, 0.5);
    }

    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        // Stimulate the apex
        if (mpMesh->GetNode(node)->rGetLocation()[0] > 0.94)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver,mpStimulus, mpZeroStimulus);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver,mpZeroStimulus, mpZeroStimulus);
        }
    }

    ~PointStimulusHeartCellFactory(void)
    {
        delete mpStimulus;
    }
};



class TestBidomainHeart : public CxxTest::TestSuite
{

public:

    void TestBidomainDg0Heart() throw (Exception)
    {
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.0025, 0.005, 0.1);                
        HeartConfig::Instance()->SetSimulationDuration(100.0);  //ms
        HeartConfig::Instance()->SetUseRelativeTolerance(5e-5);
        HeartConfig::Instance()->SetMeshFileName("heart/test/data/halfheart");
        HeartConfig::Instance()->SetOutputDirectory("BiDg0Heart");
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_Heart");
                
        PointStimulusHeartCellFactory cell_factory;
        BidomainProblem<3> bidomain_problem(&cell_factory);

        PetscOptionsSetValue("-options_table", "");

        bidomain_problem.SetWriteInfo();

        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        EventHandler::Headings();
        EventHandler::Report();
    }

    // Creates data for the following test
    void TestPermuteWithMetisBinaries()
    {
        EXIT_IF_SEQUENTIAL;
        
        TrianglesMeshReader<3,3> mesh_reader("heart/test/data/halfheart");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        unsigned num_procs = PetscTools::NumProcs();
        mesh.PermuteNodesWithMetisBinaries(num_procs);

        TrianglesMeshWriter<3,3> mesh_writer("","halfheart_metis");
        mesh_writer.WriteFilesUsingMesh(mesh);
    }


    void TestBidomainDg0HeartMetis() throw (Exception)
    {
        EXIT_IF_SEQUENTIAL;
        
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75, 1.75, 1.75));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(7.0, 7.0, 7.0));                

        HeartConfig::Instance()->SetPrintingTimeStep(0.1);        
        HeartConfig::Instance()->SetPdeTimeStep(0.005);       
        HeartConfig::Instance()->SetOdeTimeStep(0.0025);
        HeartConfig::Instance()->SetSimulationDuration(100.0);  //ms
        HeartConfig::Instance()->SetUseRelativeTolerance(5e-5);
        HeartConfig::Instance()->SetOutputDirectory("BiDg0HeartMetis");
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_HeartMetis");

        PointStimulusHeartCellFactory cell_factory;
        BidomainProblem<3> bidomain_problem(&cell_factory);

        // Data created in TestPermuteWithMetisBinaries
        OutputFileHandler handler("");
        std::string metis_mesh = handler.GetOutputDirectoryFullPath("") + "halfheart_metis";
        std::string nodes_file = handler.GetOutputDirectoryFullPath("") + "metis.mesh.nodesperproc";

        HeartConfig::Instance()->SetMeshFileName(metis_mesh);//"heart/test/data/halfheart_metis");
        bidomain_problem.SetNodesPerProcessorFilename(nodes_file);

        //PetscOptionsSetValue("-ksp_type", "symmlq");
        //PetscOptionsSetValue("-pc_type", "bjacobi");
        //PetscOptionsSetValue("-log_summary", "");
        //PetscOptionsSetValue("-ksp_monitor", "");
        PetscOptionsSetValue("-options_table", "");

        bidomain_problem.SetWriteInfo();

        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        EventHandler::Headings();
        EventHandler::Report();
    }
};

#endif //_TESTBIDOMAINHEART_HPP_
