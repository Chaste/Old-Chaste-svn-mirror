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

#include "CardiacSimulationArchiver.hpp" // Needs to be before other Chaste code

#include "Exception.hpp"
#include "DistributedVector.hpp"
#include "DistributedVectorFactory.hpp"
#include "ArchiveOpener.hpp"
#include "AbstractCardiacCell.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "BackwardEulerFoxModel2002Modified.hpp"
#include "BidomainProblem.hpp"
#include "MonodomainProblem.hpp"
#include "CompareHdf5ResultsFiles.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestCardiacSimulationArchiver : public CxxTest::TestSuite
{
private:

    std::vector<double> mSolutionReplicated1d2ms;///<Used to test differences between tests

public:

    /*
     *  Simple bidomain simulation to test against in TestArchivingWithHelperClass below
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
        
        {
            // Coverage of "couldn't find file" exception
            BidomainProblem<1> *p_bidomain_problem;            
            TS_ASSERT_THROWS_CONTAINS(p_bidomain_problem = CardiacSimulationArchiver<BidomainProblem<1> >::Load("missing_directory"),
                                      "Cannot load main archive file:");
        }

    }
    
    /**
     *  Test used to generate data for the acceptance test resume_bidomain. We run the same simulation as in save_bidomain
     *  and archive it. resume_bidomain will load it and resume the simulation.
     * 
     *  If the archiving format changes, both acceptance tests (save_bidomain and resume_bidomain) and the second part of this 
     *  test will fail. Do the following to fix them.
     * 
     *  Updated save_bidomain/ChasteResults_10ms_arch_0.chaste
     *    This can be easily done with texttest GUI. Run save_bidomain test to find that ChasteResults_10ms_arch_0.chaste
     *    contains differences (but NO other file). Highlight it and use the Save button on the window top left part.    
     * 
     *  Change into the output directory
     *    cd /tmp/chaste/testoutput
     * 
     *  Copy the archive generated.
     *    cp -r save_bidomain/ ~/eclipse/workspace/Chaste/apps/texttest/chaste/resume_bidomain/
     * 
     *  Copy the h5 results directory so it can be extended.
     *    cp -r SaveBidomain/ ~/eclipse/workspace/Chaste/apps/texttest/chaste/resume_bidomain/
     */      
    void TestGenerateResultsForResumeBidomain()
    {
        EXIT_IF_PARALLEL;
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
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_3d"); ///\todo it's not LR91

        // This cell factory should apply the same stimulus described in the xml config file.
        PlaneStimulusCellFactory<BackwardEulerFoxModel2002Modified, 3> cell_factory(-80000.0, 1.0);
        BidomainProblem<3> bidomain_problem( &cell_factory );

        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        CardiacSimulationArchiver<BidomainProblem<3> >::Save(bidomain_problem, "save_bidomain", false);
        
        {
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(
                "apps/texttest/chaste/resume_bidomain/save_bidomain/",
                "save_bidomain.arch",
                false);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();
            
            BidomainProblem<3> *p_bidomain_problem;
            TS_ASSERT_THROWS_NOTHING((*p_arch) >> p_bidomain_problem);
            //Double check that the acceptance test archive really is valid
            TS_ASSERT_EQUALS(p_bidomain_problem->mMeshFilename, "");
            TS_ASSERT_EQUALS(p_bidomain_problem->mPrintOutput, true);
            TS_ASSERT_EQUALS(p_bidomain_problem->mNodesToOutput.size(), 0u);
            TS_ASSERT_EQUALS(p_bidomain_problem->mCurrentTime, 10.0);
            TS_ASSERT_EQUALS(p_bidomain_problem->mArchiveKSP, false);
            delete p_bidomain_problem;
        }
    }    


    /**
     *  Test used to generate data for the acceptance test resume_monodomain. We run the same simulation as in save_monodomain
     *  and archive it. resume_monodomain will load it and resume the simulation.
     * 
     *  If the archiving format changes, both acceptance tests (save_monodomain and resume_monodomain) and the second part of this 
     *  test will fail. Do the following to fix them.
     * 
     *  Updated save_monodomain/ChasteResults_10ms_arch_0.chaste
     *    This can be easily done with texttest GUI. Run save_monodomain test to find that ChasteResults_10ms_arch_0.chaste
     *    contains differences (but NO other file). Highlight it and use the Save button on the window top left part.    
     * 
     *  Change into the output directory
     *    cd /tmp/chaste/testoutput
     * 
     *  Copy the archive generated.
     *    cp -r save_monodomain/ ~/eclipse/workspace/Chaste/apps/texttest/chaste/resume_monodomain/
     * 
     *  Copy the h5 results directory so it can be extended.
     *    cp -r SaveMonodomain/ ~/eclipse/workspace/Chaste/apps/texttest/chaste/resume_monodomain/
     */      
    void TestGenerateResultsForResumeMonodomain()
    {
        EXIT_IF_PARALLEL;
        HeartConfig::Instance()->SetUseFixedSchemaLocation(true);
        HeartConfig::Instance()->SetParametersFile("apps/texttest/chaste/save_monodomain/ChasteParameters.xml");

        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetSimulationDuration(), 10.0);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetDefaultIonicModel(),
                         cp::ionic_models_available_type::Fox2002BackwardEuler);

        HeartConfig::Instance()->SetOutputDirectory("SaveMonodomain");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_2d"); ///\todo it's not LR91

        // This cell factory should apply the same stimulus described in the xml config file.
        PlaneStimulusCellFactory<BackwardEulerFoxModel2002Modified, 2> cell_factory(-600000.0, 1.0);
        MonodomainProblem<2> monodomain_problem( &cell_factory );

        monodomain_problem.Initialise();
        monodomain_problem.Solve();

        CardiacSimulationArchiver<MonodomainProblem<2> >::Save(monodomain_problem, "save_monodomain", false);
        
        {
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(
                "apps/texttest/chaste/resume_monodomain/save_monodomain/",
                "save_monodomain.arch",
                false);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();
            
            MonodomainProblem<2> *p_monodomain_problem;
            TS_ASSERT_THROWS_NOTHING((*p_arch) >> p_monodomain_problem);
            //Double check that the acceptance test archive really is valid
            TS_ASSERT_EQUALS(p_monodomain_problem->mMeshFilename, "");
            TS_ASSERT_EQUALS(p_monodomain_problem->mPrintOutput, true);
            TS_ASSERT_EQUALS(p_monodomain_problem->mNodesToOutput.size(), 0u);
            TS_ASSERT_EQUALS(p_monodomain_problem->mCurrentTime, 10.0);
            TS_ASSERT_EQUALS(p_monodomain_problem->mArchiveKSP, false);
            delete p_monodomain_problem;
        }
    }

    /**
     * Run this in parallel (build=_3) to create the archive for TestLoadAsSequential.
     * Then do
     *   cd /tmp/chaste/testoutput/TestCreateArchiveForLoadAsSequential
     *   TO_DIR="$HOME/eclipse/workspace/Chaste/heart/test/data/checkpoint_migration/"
     *   for f in T*; do cp $f $TO_DIR/Test${f:20:${#f}}; done
     *   for f in [^T]*; do cp $f $TO_DIR/$f; done
     * 
     * Sets up a simulation and archives it without solving at all.
     * 
     * When running sequentially, this creates an archive we can compare with
     * that produced by the next test.
     * 
     * Generates a 3d cube mesh with 125 nodes, corners at (0,0,0) and (1,1,1)
     * with nodal spacing of 0.25cm.
     */
    void TestCreateArchiveForLoadAsSequential() throw (Exception)
    {
        std::string directory = "TestCreateArchiveForLoadAsSequential";
        HeartConfig::Instance()->Reset();
        HeartConfig::Instance()->SetSlabDimensions(1, 1, 1, 0.25);
        HeartConfig::Instance()->SetSimulationDuration(0.2);
        HeartConfig::Instance()->SetOutputDirectory(directory);
        HeartConfig::Instance()->SetOutputFilenamePrefix("simulation");
        // We want the numbers matching in any h5 files:
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-6);
        
        PlaneStimulusCellFactory<BackwardEulerFoxModel2002Modified, 3> cell_factory(-80000.0, 1.0);
        BidomainProblem<3> bidomain_problem( &cell_factory );
        
        bidomain_problem.Initialise();
        
        CardiacSimulationArchiver<BidomainProblem<3> >::Save(bidomain_problem, directory);
    }

    /**
     * #1159 - the first part of migrating a checkpoint to a different number of processes.
     */
    void TestLoadAsSequential() throw (Exception)
    {
        // We can only load simulations from CHASTE_TEST_OUTPUT, so copy the archives there
        std::string source_directory = "heart/test/data/checkpoint_migration/";
        std::string archive_directory = "TestLoadAsSequential";
        OutputFileHandler handler(archive_directory); // Clear the target directory
        if (PetscTools::AmMaster())
        {
            EXPECT0(system, "cp " + source_directory + "* " + handler.GetOutputDirectoryFullPath());
        }
        
        if (PetscTools::IsSequential())
        {
            // Do the migration to sequential
            BidomainProblem<3>* p_problem = CardiacSimulationArchiver<BidomainProblem<3> >::LoadAsSequential(archive_directory);
            
            // Test we have the right data
            const unsigned num_cells = 125u;
            TS_ASSERT_EQUALS(p_problem->mMeshFilename, "");
            TS_ASSERT_EQUALS(p_problem->mPrintOutput, true);
            TS_ASSERT_EQUALS(p_problem->mNodesToOutput.size(), 0u);
            TS_ASSERT_EQUALS(p_problem->mCurrentTime, 0.0);
            TS_ASSERT_EQUALS(p_problem->GetPde()->GetCellsDistributed().size(), num_cells);
            TS_ASSERT_EQUALS(p_problem->rGetMesh().GetNumAllNodes(), num_cells);
            TS_ASSERT_EQUALS(p_problem->rGetMesh().GetNumNodes(), num_cells);
            TS_ASSERT_EQUALS(&(p_problem->rGetMesh()), p_problem->GetPde()->pGetMesh());
            
            // All cells at x=0 should have a SimpleStimulus(-80000, 1).
            // Also all cells should be at initial conditions.
            std::vector<double> inits = p_problem->GetPde()->GetCardiacCell(0)->GetInitialConditions();
            for (unsigned i=0; i<num_cells; i++)
            {
                AbstractCardiacCell* p_cell = p_problem->GetPde()->GetCardiacCell(i);
                const boost::shared_ptr<AbstractStimulusFunction> p_stim = p_cell->GetStimulusFunction();
                double x = p_problem->rGetMesh().GetNode(i)->GetPoint()[0];
                
                if (x*x < 1e-10)
                {
                    // Stim exists
                    TS_ASSERT_DELTA(p_cell->GetStimulus(0.0), -80000.0, 1e-10);
                    TS_ASSERT_DELTA(p_cell->GetStimulus(1.0), -80000.0, 1e-10);
                    TS_ASSERT_DELTA(p_cell->GetStimulus(1.001), 0.0, 1e-10);
                    TS_ASSERT_DELTA(p_cell->GetStimulus(-1e-10), 0.0, 1e-10);
                }
                else
                {
                    // No stim
                    TS_ASSERT_DELTA(p_cell->GetStimulus(0.0), 0.0, 1e-10);
                    TS_ASSERT_DELTA(p_cell->GetStimulus(1.0), 0.0, 1e-10);
                    TS_ASSERT_DELTA(p_cell->GetStimulus(1.001), 0.0, 1e-10);
                    TS_ASSERT_DELTA(p_cell->GetStimulus(-1e-10), 0.0, 1e-10);
                }
                
                // All cells should be at initial conditions
                std::vector<double>& r_state = p_cell->rGetStateVariables();
                TS_ASSERT_EQUALS(r_state.size(), inits.size());
                for (unsigned j=0; j<r_state.size(); j++)
                {
                    TS_ASSERT_DELTA(r_state[j], inits[j], 1e-10);
                }
            }
            
            // Test bccs - none defined in this problem
            TS_ASSERT(! p_problem->mpDefaultBoundaryConditionsContainer);
            TS_ASSERT(! p_problem->mpBoundaryConditionsContainer);
            
            // Save it to a sequential archive
            CardiacSimulationArchiver<BidomainProblem<3> >::Save(*p_problem, archive_directory);
            
            // Compare with the archive from the previous test
            std::string ref_archive_dir = "TestCreateArchiveForLoadAsSequential";
            std::string ref_archive = handler.GetChasteTestOutputDirectory() + ref_archive_dir + "/" + ref_archive_dir + ".arch";
            std::string my_archive = handler.GetOutputDirectoryFullPath() + archive_directory + ".arch";
            EXPECT0(system, "diff " + ref_archive + " " + my_archive);
            // This will differ because we get extra copies of the ODE solver and intracellular stimulus objects
            // (one per original process which had them).
            //EXPECT0(system, "diff " + ref_archive + ".0 " + my_archive + ".0");
            EXPECT0(system, "diff -I 'serialization::archive' " + source_directory + "reference_0_archive " + my_archive + ".0");
            
            // Simulate this problem
            // Change output directory to avoid over-writing original results & archives
            HeartConfig::Instance()->SetOutputDirectory(ref_archive_dir + "/mig1");
            p_problem->Solve();
            // Copy the results vector
            Vec migrated_soln = p_problem->GetSolution();
            Vec migrated_soln_copy;
            VecDuplicate(migrated_soln, &migrated_soln_copy);
            VecCopy(migrated_soln, migrated_soln_copy);
            // and destroy the problem, so we don't get confusion from 2 problems at the same time
            delete p_problem;
            
            // Compare the results with simulating the archive from the previous test
            BidomainProblem<3>* p_orig_problem = CardiacSimulationArchiver<BidomainProblem<3> >::Load(ref_archive_dir);
            HeartConfig::Instance()->SetOutputDirectory(ref_archive_dir + "/orig");
            p_orig_problem->Solve();
            DistributedVector orig_soln = p_orig_problem->GetSolutionDistributedVector();
            DistributedVector::Stripe orig_vm(orig_soln, 0);
            DistributedVector::Stripe orig_phie(orig_soln, 1);
            DistributedVector migrated_soln_1(migrated_soln_copy, orig_soln.GetFactory());
            DistributedVector::Stripe migrated_vm_1(migrated_soln_1, 0);
            DistributedVector::Stripe migrated_phie_1(migrated_soln_1, 1);
            for (DistributedVector::Iterator index = migrated_soln_1.Begin();
                 index != migrated_soln_1.End();
                 ++index)
            {
                TS_ASSERT_DELTA(orig_vm[index], migrated_vm_1[index], 1e-8);
                TS_ASSERT_DELTA(orig_phie[index], migrated_phie_1[index], 1e-8);
            }
            delete p_orig_problem;
            
            // Now try loading the migrated simulation that we saved above
            p_problem = CardiacSimulationArchiver<BidomainProblem<3> >::Load(archive_directory);
            HeartConfig::Instance()->SetOutputDirectory(ref_archive_dir + "/mig2");
            p_problem->Solve();
            // and again compare the results
            DistributedVector migrated_soln_2 = p_problem->GetSolutionDistributedVector();
            DistributedVector::Stripe migrated_vm_2(migrated_soln_2, 0);
            DistributedVector::Stripe migrated_phie_2(migrated_soln_2, 1);
            for (DistributedVector::Iterator index = migrated_soln_1.Begin();
                 index != migrated_soln_1.End();
                 ++index)
            {
                TS_ASSERT_DELTA(migrated_vm_2[index], migrated_vm_1[index], 1e-8);
                TS_ASSERT_DELTA(migrated_phie_2[index], migrated_phie_1[index], 1e-8);
            }
            delete p_problem;
            
            // Cover exception
            TS_ASSERT_THROWS_CONTAINS(p_problem = CardiacSimulationArchiver<BidomainProblem<3> >::LoadAsSequential("non_existent_dir"),
                                      "Cannot load main archive file: ");
        }
        else
        {
            BidomainProblem<3>* p_problem;
            TS_ASSERT_THROWS_THIS(p_problem = CardiacSimulationArchiver<BidomainProblem<3> >::LoadAsSequential(archive_directory),
                                  "Cannot load sequentially when running in parallel.");
        }
    }
};

#endif /*TESTCARDIACSIMULATIONARCHIVER_HPP_*/
