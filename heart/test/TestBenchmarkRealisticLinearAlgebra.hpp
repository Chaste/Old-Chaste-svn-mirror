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

#ifndef TESTREALISTICLINEARALGEBRA_HPP_
#define TESTREALISTICLINEARALGEBRA_HPP_

#include <cxxtest/TestSuite.h>

#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "BidomainProblem.hpp"
#include "PetscSetupAndFinalize.hpp"

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
    
    
    
class TestBenchmarkRealisticLinearAlgebra : public CxxTest::TestSuite
{
private:
    void SetParameters()
    {
        // Timesteps and simulation duration
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.0025, 0.005, 0.1);                
        HeartConfig::Instance()->SetSimulationDuration(1.0);  //ms

        // Input cardiac mesh
        HeartConfig::Instance()->SetMeshFileName("heart/test/data/halfheart");
    
        // Output directory
        HeartConfig::Instance()->SetOutputDirectory("BenchmarkRealisticLA");
        
        // Solver and preconditioner selection through Chaste parameter system.
        // Not all the possible methods can be selected via HeartConfig. If an error like the following
        // happens, see next comment.
        //    heart/test/TestRealisticLinearAlgebra.hpp:105: Error: Test failed: 
        //    Chaste error: heart/src/problem/HeartConfig.cpp:866: Unknown solver type provided        
        HeartConfig::Instance()->SetKSPSolver("symmlq");
        HeartConfig::Instance()->SetKSPPreconditioner("bjacobi");

        // In the case you want to select a solver or preconditioner not supported in HeartConfig,
        // you should talk to the KSP object directly. Uncomment and modify accordingly
//        PetscOptionsSetValue("-ksp_type", "bicg");
//        PetscOptionsSetValue("-pc_type", "asm");
                
        // Traces KSP solution (# of iterations, residual, etc)
        PetscOptionsSetValue("-ksp_monitor", "");
        
        // Enables extra logging (# of flops, messages, reductions, etc)
//        PetscOptionsSetValue("-log_summary", "");        
    }
    
public:

    void TestBenchmarkBidomain() throw (Exception)
    {
        SetParameters();

        // Output filename for this particular test
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLA");
        
        PointStimulusHeartCellFactory cell_factory;
        BidomainProblem<3> bidomain_problem(&cell_factory);
    
        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        EventHandler::Headings();
        EventHandler::Report();
    }
};

#endif /*TESTREALISTICLINEARALGEBRA_HPP_*/
