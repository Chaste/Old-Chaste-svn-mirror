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
    PointStimulusHeartCellFactory(double timeStep) : AbstractCardiacCellFactory<3>(timeStep)
    {
        mpStimulus = new SimpleStimulus(-1000.0*100, 0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {      
        // Stimulate the apex
        if (mpMesh->GetNode(node)->rGetLocation()[0] > 0.94)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpStimulus, mpZeroStimulus);            
        }        
        else
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpZeroStimulus, mpZeroStimulus);
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
        double pde_time_step = 0.005;  // ms
        double ode_time_step = 0.0025; // ms
        double end_time = 100;        // ms
        double printing_time_step = 0.1;
        
        PointStimulusHeartCellFactory cell_factory(ode_time_step);
        BidomainProblem<3> bidomain_problem(&cell_factory);
        
        bidomain_problem.SetMeshFilename("heart/test/data/halfheart");
        bidomain_problem.SetOutputDirectory("BiDg0Heart");
        bidomain_problem.SetOutputFilenamePrefix("BidomainLR91_Heart");
        
        bidomain_problem.SetEndTime(end_time);
        bidomain_problem.SetPdeAndPrintingTimeSteps(pde_time_step, printing_time_step);
        
        //bidomain_problem.SetLinearSolverRelativeTolerance(5e-7);
        bidomain_problem.SetLinearSolverAbsoluteTolerance(5e-3);        
        //PetscOptionsSetValue("-ksp_type", "symmlq");
        //PetscOptionsSetValue("-pc_type", "bjacobi");
        //PetscOptionsSetValue("-log_summary", "");
        //PetscOptionsSetValue("-ksp_monitor", "");
        PetscOptionsSetValue("-options_table", "");
        
        bidomain_problem.SetWriteInfo();

        bidomain_problem.SetIntracellularConductivities(Create_c_vector(1.75, 1.75, 1.75));
        bidomain_problem.SetExtracellularConductivities(Create_c_vector(7.0, 7.0, 7.0));
        
        bidomain_problem.Initialise();
        bidomain_problem.Solve();
        
        EventHandler::Headings();
        EventHandler::Report();
    }

    // Creates data for the following test
    void TestPermuteWithMetisBinaries()
    {
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
        double pde_time_step = 0.005;  // ms
        double ode_time_step = 0.0025; // ms
        double end_time = 100;        // ms
        double printing_time_step = 0.1;
        
        PointStimulusHeartCellFactory cell_factory(ode_time_step);
        BidomainProblem<3> bidomain_problem(&cell_factory);
        
        // Data created in TestPermuteWithMetisBinaries
        OutputFileHandler handler("");
        std::string metis_mesh = handler.GetOutputDirectoryFullPath("") + "halfheart_metis";
        std::string nodes_file = handler.GetOutputDirectoryFullPath("") + "metis.mesh.nodesperproc";
        
        bidomain_problem.SetMeshFilename(metis_mesh);//"heart/test/data/halfheart_metis");
        bidomain_problem.SetNodesPerProcessorFilename(nodes_file);
        bidomain_problem.SetOutputDirectory("BiDg0HeartMetis");
        bidomain_problem.SetOutputFilenamePrefix("BidomainLR91_HeartMetis");
        
        bidomain_problem.SetEndTime(end_time);
        bidomain_problem.SetPdeAndPrintingTimeSteps(pde_time_step, printing_time_step);
        
        //bidomain_problem.SetLinearSolverRelativeTolerance(5e-7);
        bidomain_problem.SetLinearSolverAbsoluteTolerance(5e-3);
        //PetscOptionsSetValue("-ksp_type", "symmlq");
        //PetscOptionsSetValue("-pc_type", "bjacobi");
        //PetscOptionsSetValue("-log_summary", "");
        //PetscOptionsSetValue("-ksp_monitor", "");
        PetscOptionsSetValue("-options_table", "");
        
        bidomain_problem.SetWriteInfo();

        bidomain_problem.SetIntracellularConductivities(Create_c_vector(1.75, 1.75, 1.75));
        bidomain_problem.SetExtracellularConductivities(Create_c_vector(7.0, 7.0, 7.0));        
        
        bidomain_problem.Initialise();
        bidomain_problem.Solve();
                
        EventHandler::Headings();
        EventHandler::Report();        
    }
};

#endif //_TESTBIDOMAINHEART_HPP_
