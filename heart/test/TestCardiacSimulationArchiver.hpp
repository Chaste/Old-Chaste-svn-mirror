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

#ifndef TESTCARDIACSIMULATIONARCHIVER_HPP_
#define TESTCARDIACSIMULATIONARCHIVER_HPP_

#include <cxxtest/TestSuite.h>

#include "CardiacSimulationArchiver.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "BidomainProblem.hpp"
#include "MonodomainProblem.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "ParallelTetrahedralMesh.hpp"
#include "CompareHdf5ResultsFiles.hpp"
#include "BackwardEulerFoxModel2002Modified.hpp"

class TestCardiacSimulationArchiver : public CxxTest::TestSuite
{
private:

    std::vector<double> mSolutionReplicated1d2ms;///<Used to test differences between tests

public:

    /*
     *  Simple bidomain simulation to test against in the archiving tests below
     */
    void TestSimpleBidomain1D() throw(Exception)
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);        
        
        HeartConfig::Instance()->SetSimulationDuration(2.0); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("BidomainSimple1d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_1d");

        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1> cell_factory;
        BidomainProblem<1> bidomain_problem( &cell_factory );

        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        // check some voltages
        ReplicatableVector solution_replicated(bidomain_problem.GetSolution());

        double atol=5e-3;

        TS_ASSERT_DELTA(solution_replicated[1], -16.4861, atol);
        TS_ASSERT_DELTA(solution_replicated[2], 22.8117, atol);
        TS_ASSERT_DELTA(solution_replicated[3], -16.4893, atol);
        TS_ASSERT_DELTA(solution_replicated[5], -16.5617, atol);
        TS_ASSERT_DELTA(solution_replicated[7], -16.6761, atol);
        TS_ASSERT_DELTA(solution_replicated[9], -16.8344, atol);
        TS_ASSERT_DELTA(solution_replicated[10], 25.3148, atol);
        
        for (unsigned index=0; index<solution_replicated.GetSize(); index++)
        {
            mSolutionReplicated1d2ms.push_back(solution_replicated[index]);
        }

    }

    /*
     *  Same as TestArchiving but through the helper class
     */    
    void TestArchivingWithHelperClass()
    {
        // Save
        {
            HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
            HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.0005));
            //HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
            HeartConfig::Instance()->SetOutputDirectory("BiProblemArchiveHelper");
            HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_1d");
            HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
            HeartConfig::Instance()->SetCapacitance(1.0);
    
            PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1> cell_factory;
            BidomainProblem<1> bidomain_problem( &cell_factory );

            /// \todo: Make this test pass if the mesh is set via HeartConfig
            TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1mm_10_elements");
            ParallelTetrahedralMesh<1,1> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);
            bidomain_problem.SetMesh(&mesh);
    
            bidomain_problem.Initialise();
            HeartConfig::Instance()->SetSimulationDuration(1.0); //ms
            bidomain_problem.Solve();
                        
            CardiacSimulationArchiver<BidomainProblem<1> >::Save(bidomain_problem, "bidomain_problem_archive_helper", false);
        }

        // Load
        {
            BidomainProblem<1> *p_bidomain_problem;            
            p_bidomain_problem = CardiacSimulationArchiver<BidomainProblem<1> >::Load("bidomain_problem_archive_helper");

            HeartConfig::Instance()->SetSimulationDuration(2.0); //ms
            p_bidomain_problem->Solve();
    
            // check some voltages
            ReplicatableVector solution_replicated(p_bidomain_problem->GetSolution());
            double atol=5e-3;
            TS_ASSERT_DELTA(solution_replicated[1], -16.4861, atol);
            TS_ASSERT_DELTA(solution_replicated[2], 22.8117, atol);
            TS_ASSERT_DELTA(solution_replicated[3], -16.4893, atol);
            TS_ASSERT_DELTA(solution_replicated[5], -16.5617, atol);
            TS_ASSERT_DELTA(solution_replicated[7], -16.6761, atol);
            TS_ASSERT_DELTA(solution_replicated[9], -16.8344, atol);
            TS_ASSERT_DELTA(solution_replicated[10], 25.3148, atol);        

            for (unsigned index=0; index<solution_replicated.GetSize(); index++)
            {
                //Shouldn't differ from the original run at all
                TS_ASSERT_DELTA(solution_replicated[index], mSolutionReplicated1d2ms[index],  5e-11);
            }
            // check output file contains results for the whole simulation
            TS_ASSERT(CompareFilesViaHdf5DataReader("BiProblemArchiveHelper", "BidomainLR91_1d", true,
                                                    "BidomainSimple1d", "BidomainLR91_1d", true));
            
            // Free memory
            delete p_bidomain_problem;            
        }

    }
    
    /*
     *  Test used to generate the results for the acceptance test save_bidomain.
     * 
     *  The idea is that the binary should generate the same results that this test.
     * 
     *  If the archiving format changes you will need to update two convergence tests based on the result of this test:
     * 
     *  save_bidomain:
     *    h5dump save_bidomain/AbstractCardiacProblem_mSolution.h5 > ~/eclipse/workspace/Chaste/apps/texttest/chaste/save_bidomain/AbstractCardiacProblem_mSolution_h5.chaste
     *    cp save_bidomain.arch.0 ~/eclipse/workspace/Chaste/apps/texttest/chaste/save_bidomain/ChasteResults_10ms_arch_0.chaste
     * 
     *  resume_bidomain:
     *    cd /tmp/chaste/testoutput
     *    cp -r SaveBidomain/ ~/eclipse/workspace/Chaste/apps/texttest/chaste/resume_bidomain/
     *    cp -r save_bidomain/ ~/eclipse/workspace/Chaste/apps/texttest/chaste/resume_bidomain/
     * 
     */      
    void TestGenerateResultsForSaveBidomain()
    {
        HeartConfig::Instance()->SetUseFixedSchemaLocation(true);
        HeartConfig::Instance()->SetParametersFile("apps/texttest/chaste/save_bidomain/ChasteParameters.xml");
        // We reset the mesh filename to include the relative path
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/cube_1626_elements");

        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetSimulationDuration(), 10.0);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetMeshName(),
                         "mesh/test/data/cube_1626_elements");
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetDefaultIonicModel(),
                         cp::ionic_models_available_type::Fox2002BackwardEuler);

        HeartConfig::Instance()->SetOutputDirectory("SaveBidomain");
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_1d");      

        // This cell factory should apply the same stimulus described in the xml config file.
        PlaneStimulusCellFactory<BackwardEulerFoxModel2002Modified, 3> cell_factory(-80000.0, 1.0);
        BidomainProblem<3> bidomain_problem( &cell_factory );

        bidomain_problem.ConvertOutputToMeshalyzerFormat(true);

        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        CardiacSimulationArchiver<BidomainProblem<3> >::Save(bidomain_problem, "save_bidomain", false);
        
        /// \todo: add a test that fails if the archiving format is modified.
    }
    
    /*
     *  Test used to generate the results for the acceptance test resume_bidomain.
     * 
     *  Here we run the whole simulation (20ms) and check
     *
     *  If the archiving format changes you will need to update two convergence tests based on the result of this test:
     * 
     *  resume_bidomain:
     *    cd /tmp/chaste/testoutput/ResumeBidomain
     *    h5dump BidomainLR91_1d.h5 > ~/eclipse/workspace/Chaste/apps/texttest/chaste/resume_bidomain/results_h5.chaste
     */
    void TestGenerateResultsForResumeBidomain()
    {
        HeartConfig::Instance()->SetParametersFile("apps/texttest/chaste/save_bidomain/ChasteParameters.xml");
        // We reset the mesh filename to include the relative path
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/cube_1626_elements");

        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetSimulationDuration(), 10.0);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetMeshName(),
                         "mesh/test/data/cube_1626_elements");
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetDefaultIonicModel(),
                         cp::ionic_models_available_type::Fox2002BackwardEuler);

        // We simulate for 20ms to compare the output
        HeartConfig::Instance()->SetSimulationDuration(20.0);

        HeartConfig::Instance()->SetOutputDirectory("ResumeBidomain");
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_1d");      

        // This cell factory should apply the same stimulus described in the xml config file.
        PlaneStimulusCellFactory<BackwardEulerFoxModel2002Modified, 3> cell_factory(-80000.0, 1.0);
        BidomainProblem<3> bidomain_problem( &cell_factory );

        bidomain_problem.ConvertOutputToMeshalyzerFormat(true);

        bidomain_problem.Initialise();
        bidomain_problem.Solve();
        
        /// \todo: add a test that fails if the archiving format is modified.
    }


    /*
     *  Test used to generate the results for the acceptance test save_monodomain.
     * 
     *  The idea is that the binary should generate the same results that this test.
     * 
     *  If the archiving format changes you will need to update two convergence tests based on the result of this test:
     * 
     *  save_monodomain:
     *    cd /tmp/chaste/testoutput
     *    h5dump save_monodomain/AbstractCardiacProblem_mSolution.h5 > ~/eclipse/workspace/Chaste/apps/texttest/chaste/save_monodomain/AbstractCardiacProblem_mSolution_h5.chaste
     *    cp save_monodomain/save_monodomain.arch.0 ~/eclipse/workspace/Chaste/apps/texttest/chaste/save_monodomain/ChasteResults_10ms_arch_0.chaste
     * 
     *  resume_monodomain:
     *    cd /tmp/chaste/testoutput
     *    cp -r SaveMonodomain/ ~/eclipse/workspace/Chaste/apps/texttest/chaste/resume_monodomain/
     *    cp -r save_monodomain/ ~/eclipse/workspace/Chaste/apps/texttest/chaste/resume_monodomain/
     * 
     */      
    void TestGenerateResultsForSaveMonodomain()
    {
        HeartConfig::Instance()->SetUseFixedSchemaLocation(true);
        HeartConfig::Instance()->SetParametersFile("apps/texttest/chaste/save_monodomain/ChasteParameters.xml");
        // We reset the mesh filename to include the relative path
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/cube_1626_elements");

        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetSimulationDuration(), 10.0);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetDefaultIonicModel(),
                         cp::ionic_models_available_type::Fox2002BackwardEuler);

        HeartConfig::Instance()->SetOutputDirectory("SaveMonodomain");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_1d");      

        // This cell factory should apply the same stimulus described in the xml config file.
        PlaneStimulusCellFactory<BackwardEulerFoxModel2002Modified, 2> cell_factory(-600000.0, 1.0);
        MonodomainProblem<2> monodomain_problem( &cell_factory );

        monodomain_problem.ConvertOutputToMeshalyzerFormat(true);

        monodomain_problem.Initialise();
        monodomain_problem.Solve();

        CardiacSimulationArchiver<MonodomainProblem<2> >::Save(monodomain_problem, "save_monodomain", false);
        
        /// \todo: add a test that fails if the archiving format is modified.
    }
    
    /*
     *  Test used to generate the results for the acceptance test resume_monodomain.
     * 
     *  Here we run the whole simulation (20ms) and check
     *
     *  If the archiving format changes you will need to update two convergence tests based on the result of this test:
     * 
     *  resume_monodomain:
     *    cd /tmp/chaste/testoutput/ResumeMonodomain
     *    h5dump MonodomainLR91_1d.h5 > ~/eclipse/workspace/Chaste/apps/texttest/chaste/resume_monodomain/resumed_mono_h5.chaste
     */
    void TestGenerateResultsForResumeMonodomain()
    {
        HeartConfig::Instance()->SetParametersFile("apps/texttest/chaste/save_monodomain/ChasteParameters.xml");
        // We reset the mesh filename to include the relative path
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/cube_1626_elements");

        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetSimulationDuration(), 10.0);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetDefaultIonicModel(),
                         cp::ionic_models_available_type::Fox2002BackwardEuler);

        // We simulate for 20ms to compare the output
        HeartConfig::Instance()->SetSimulationDuration(20.0);

        HeartConfig::Instance()->SetOutputDirectory("ResumeMonodomain");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_1d");      

        // This cell factory should apply the same stimulus described in the xml config file.
        PlaneStimulusCellFactory<BackwardEulerFoxModel2002Modified, 2> cell_factory(-600000.0, 1.0);
        MonodomainProblem<2> monodomain_problem( &cell_factory );

        monodomain_problem.ConvertOutputToMeshalyzerFormat(true);

        monodomain_problem.Initialise();
        monodomain_problem.Solve();
        
        /// \todo: add a test that fails if the archiving format is modified.
    }

};

#endif /*TESTCARDIACSIMULATIONARCHIVER_HPP_*/
