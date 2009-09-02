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


#ifndef _TESTMONODOMAINPROBLEM_HPP_
#define _TESTMONODOMAINPROBLEM_HPP_


#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <vector>
#include "MonodomainProblem.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "FaberRudy2000Version3.hpp"
#include "Hdf5DataReader.hpp"
#include "ReplicatableVector.hpp"
#include "CheckMonoLr91Vars.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "TetrahedralMesh.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"


/*
 *  This cell factory introduces a stimulus in the very centre of mesh/test/data/2D_0_to_1mm_400_elements.
 *  This is node 60 at (0.5, 0.5)
 */
class PointStimulus2dCellFactory : public AbstractCardiacCellFactory<2>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;
    unsigned mFoundMiddlePoint;


public:
    PointStimulus2dCellFactory()
        : AbstractCardiacCellFactory<2>(),
          mpStimulus(new SimpleStimulus(-6000.0, 0.5)),
          mFoundMiddlePoint(0)
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        ChastePoint<2> location = GetMesh()->GetNode(node)->GetPoint();

        if (fabs(location[0]-0.05)<1e-6 && fabs(location[1]-0.05)<1e-6)
        {
            mFoundMiddlePoint++;
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpStimulus);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpZeroStimulus);
        }
    }

    void FinaliseCellCreation(std::vector<AbstractCardiacCell* >* pCellsDistributed, unsigned lo, unsigned hi)
    {
        unsigned found_middle_point_reduced;
        MPI_Allreduce(&mFoundMiddlePoint, &found_middle_point_reduced, 1, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);
        assert(found_middle_point_reduced == 1); // Only 1 cell should be stimulated
    }
};



class TestMonodomainProblem : public CxxTest::TestSuite
{
public:
    void tearDown()
    {
        HeartConfig::Reset();
    }

    // Solve on a 1D string of cells, 1mm long with a space step of 0.1mm.
    void TestMonodomainProblem1D() throw(Exception)
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetSimulationDuration(2.0); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("MonoProblem1d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_1d");

        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1> cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );

        monodomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        monodomain_problem.Solve();

        // test whether voltages and gating variables are in correct ranges
        CheckMonoLr91Vars<1>(monodomain_problem);

        // check some voltages
        ReplicatableVector voltage_replicated(monodomain_problem.GetSolution());
        double atol=5e-3;

        TS_ASSERT_DELTA(voltage_replicated[1], 20.7710232, atol);
        TS_ASSERT_DELTA(voltage_replicated[3], 21.5319692, atol);
        TS_ASSERT_DELTA(voltage_replicated[5], 22.9280817, atol);
        TS_ASSERT_DELTA(voltage_replicated[7], 24.0611303, atol);
        TS_ASSERT_DELTA(voltage_replicated[9], -0.770330519, atol);
        TS_ASSERT_DELTA(voltage_replicated[10], -19.2234919, atol);

        // cover get pde
        monodomain_problem.GetPde();

        // check a progress report exists
        TS_ASSERT_EQUALS(system(("ls " + OutputFileHandler::GetChasteTestOutputDirectory() + "MonoProblem1d/").c_str()), 0);
    }

    void TestMonodomainProblem1DWithRelativeTolerance() throw(Exception)
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetSimulationDuration(2.0); //ms
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);
        HeartConfig::Instance()->SetUseRelativeTolerance(1e-9);
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("MonoProblem1d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_1d");

        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1> cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );

        monodomain_problem.Initialise();

        monodomain_problem.Solve();

        // test whether voltages and gating variables are in correct ranges
        CheckMonoLr91Vars<1>(monodomain_problem);

        // check some voltages
        ReplicatableVector voltage_replicated(monodomain_problem.GetSolution());

        /// \todo: If we request "relative" tolerance we shouldn't testing in an "absolute" manner
        double atol=1e-5;
        TS_ASSERT_DELTA(voltage_replicated[1], 20.7710232, atol);
        TS_ASSERT_DELTA(voltage_replicated[3], 21.5319692, atol);
        TS_ASSERT_DELTA(voltage_replicated[5], 22.9280817, atol);
        TS_ASSERT_DELTA(voltage_replicated[7], 24.0611303, atol);
        TS_ASSERT_DELTA(voltage_replicated[9], -0.770330519, atol);
        TS_ASSERT_DELTA(voltage_replicated[10], -19.2234919, atol);

    }

    // Same as TestMonodomainProblem1D, except the 1D mesh is embedded in 3D space.
    void TestMonodomainProblem1Din3D() throw(Exception)
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetSimulationDuration(2.0); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_in_3D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("MonoProblem1din3d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_1din3d");

        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1, 3> cell_factory;
        MonodomainProblem<1,3> monodomain_problem( &cell_factory );

        monodomain_problem.ConvertOutputToMeshalyzerFormat(true);
        monodomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        monodomain_problem.Solve();

        // test whether voltages and gating variables are in correct ranges
        CheckMonoLr91Vars(monodomain_problem);

        // check some voltages
        ReplicatableVector voltage_replicated(monodomain_problem.GetSolution());
        double atol=5e-3;

        TS_ASSERT_DELTA(voltage_replicated[1], 20.7710232, atol);
        TS_ASSERT_DELTA(voltage_replicated[3], 21.5319692, atol);
        TS_ASSERT_DELTA(voltage_replicated[5], 22.9280817, atol);
        TS_ASSERT_DELTA(voltage_replicated[7], 24.0611303, atol);
        TS_ASSERT_DELTA(voltage_replicated[9], -0.770330519, atol);
        TS_ASSERT_DELTA(voltage_replicated[10], -19.2234919, atol);

        // cover get pde
        monodomain_problem.GetPde();

        // check a progress report exists
        TS_ASSERT_EQUALS(system(("ls " + OutputFileHandler::GetChasteTestOutputDirectory() + "MonoProblem1din3d/").c_str()), 0);
    }

    void TestMonodomainProblem1DWithAbsoluteTolerance() throw (Exception)
    {
        double atol = 1e-4;
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetSimulationDuration(2); //ms
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);
        HeartConfig::Instance()->SetUseAbsoluteTolerance(atol/10.0);
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("MonoProblem1d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_1d");

        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1> cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );

        monodomain_problem.Initialise();

        monodomain_problem.Solve();

        // test whether voltages and gating variables are in correct ranges
        CheckMonoLr91Vars<1>(monodomain_problem);

        // check some voltages
        ReplicatableVector voltage_replicated(monodomain_problem.GetSolution());

        TS_ASSERT_DELTA(voltage_replicated[1], 20.7710232, atol);
        TS_ASSERT_DELTA(voltage_replicated[3], 21.5319692, atol);
        TS_ASSERT_DELTA(voltage_replicated[5], 22.9280817, atol);
        TS_ASSERT_DELTA(voltage_replicated[7], 24.0611303, atol);
        TS_ASSERT_DELTA(voltage_replicated[9], -0.770330519, atol);
        TS_ASSERT_DELTA(voltage_replicated[10], -19.2234919, atol);

    }

    // Solve on a 2D 1mm by 1mm mesh (space step = 0.1mm), stimulating the left
    // edge.
    // Should behave like the 1D case, extrapolated.
    // See also TestMonodomainSlab.hpp (nightly test) for the 3D version.
    void TestMonodomainProblem2DWithEdgeStimulus() throw(Exception)
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005, 0.0005));
        HeartConfig::Instance()->SetSimulationDuration(2); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/2D_0_to_1mm_400_elements");
        HeartConfig::Instance()->SetOutputDirectory("MonoProblem2dWithEdgeStimulus");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_2dWithEdgeStimulus");

        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 2> cell_factory;

        // using the criss-cross mesh so wave propagates properly
        MonodomainProblem<2> monodomain_problem( &cell_factory );

        monodomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        monodomain_problem.Solve();

        // test whether voltages and gating variables are in correct ranges
        CheckMonoLr91Vars(monodomain_problem);

        //Since we are going to compare voltages that may be owned by
        //various processes it makes sense to replicate the data.
        Vec voltage=monodomain_problem.GetSolution();
        ReplicatableVector voltage_replicated;
        voltage_replicated.ReplicatePetscVector(voltage);
        /*
         * Test the top right node against the right one in the 1D case,
         * comparing voltage, and then test all the nodes on the right hand
         * side of the square against the top right one, comparing voltage.
         */
        bool need_initialisation = true;
        double probe_voltage=0.0;

        need_initialisation = true;

        // Test the RHS of the mesh
        AbstractTetrahedralMesh<2,2>& r_mesh=monodomain_problem.rGetMesh();
        for (AbstractTetrahedralMesh<2,2>::NodeIterator iter=r_mesh.GetNodeIteratorBegin();
             iter != r_mesh.GetNodeIteratorEnd();
             ++iter)
        {
            unsigned i=(*iter).GetIndex();
            if ((*iter).GetPoint()[0] == 0.1)
            {
                // x = 0 is where the stimulus has been applied
                // x = 0.1cm is the other end of the mesh and where we want to
                //       to test the value of the nodes

                if (need_initialisation)
                {
                    probe_voltage = voltage_replicated[i];
                    need_initialisation = false;
                }
                else
                {
                    // Tests the final voltages for all the RHS edge nodes
                    // are close to each other.
                    // This works as we are using the 'criss-cross' mesh,
                    // the voltages would vary more with a mesh with all the
                    // triangles aligned in the same direction.
                    TS_ASSERT_DELTA(voltage_replicated[i], probe_voltage, 3e-4);
                }


                // Check against 1d case - THIS TEST HAS BEEN REMOVED AS THE MESH
                // IS FINER THAN THE 1D MESH SO WE DONT EXPECT THE RESULTS TO BE THE SAME
                // TS_ASSERT_DELTA(p_voltage_array[i], -35.1363, 35*0.1);

                // test the RHS edge voltages
                // hardcoded result that looks accurate - this is a test to see
                // that nothing has changeed
                // assumes endtime = 2ms
                TS_ASSERT_DELTA(voltage_replicated[i], -59.6488, 5e-4);
            }
        }

    }


    // Solve on a 2D 1mm by 1mm mesh (space step = 0.1mm), stimulating in the
    // very centre of the mesh.
    void TestMonodomainProblem2DWithPointStimulusInTheVeryCentreOfTheMesh() throw(Exception)
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005, 0.0005));
        HeartConfig::Instance()->SetSimulationDuration(1.3); //ms - needs to be 1.3 ms to pass test
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/2D_0_to_1mm_400_elements");
        HeartConfig::Instance()->SetOutputDirectory("MonoProblem2dWithPointStimulus");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_2dWithPointStimulus");

        // To time the solve
        time_t start,end;
        double dif;
        time (&start);

        PointStimulus2dCellFactory cell_factory;

        MonodomainProblem<2> monodomain_problem( &cell_factory );

        monodomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        monodomain_problem.Solve();

        // To time the solve
        time (&end);
        dif = difftime (end,start);
        //printf ("\nSolve took %.2lf seconds. \n", dif );

        // test whether voltages and gating variables are in correct ranges
        CheckMonoLr91Vars(monodomain_problem);

        /*
         * Test that corners are 'equal', and centres of sides.
         * Irregularities in which way the triangles are oriented make
         * this rather difficult, especially since the edges are sampled
         * during the upstroke.
         */
        DistributedVector voltage = monodomain_problem.GetSolutionDistributedVector();

        // corners -> 0, 10, 110, 120
        // hardcoded result to check against
        // assumes endtime = 1.3
        unsigned corners_checked=0;
        for (DistributedVector::Iterator node_index = voltage.Begin();
             node_index!= voltage.End();
             ++node_index)
        {
            ChastePoint<2> location = monodomain_problem.rGetMesh().GetNode(node_index.Global)->GetPoint();

            if (fabs(location[0]-0.0)<1e-6 && fabs(location[1]-0.0)<1e-6) // Corner 0
            {
                TS_ASSERT_DELTA(voltage[node_index.Global], -34.3481, 1e-3);
                corners_checked++;
            }

            if (fabs(location[0]-0.1)<1e-6 && fabs(location[1]-0.0)<1e-6) // Corner 10
            {
                TS_ASSERT_DELTA(voltage[node_index.Global], -34.3481, 1e-3);
                corners_checked++;
            }

            if (fabs(location[0]-0.0)<1e-6 && fabs(location[1]-0.0)<1e-6) // Corner 110
            {
                TS_ASSERT_DELTA(voltage[node_index.Global], -34.3481, 1e-3);
                corners_checked++;
            }

            if (fabs(location[0]-0.0)<1e-6 && fabs(location[1]-0.0)<1e-6) // Corner 120
            {
                TS_ASSERT_DELTA(voltage[node_index.Global], -34.3481, 1e-3);
                corners_checked++;
            }
        }

        unsigned corners_checked_reduced;
        MPI_Allreduce(&corners_checked, &corners_checked_reduced, 1, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);
        TS_ASSERT(corners_checked_reduced==4);


        // centre of edges -> 5, 55, 65, 115
        // hardcoded result to check against
        // assumes endtime = 1.3
        unsigned edges_checked=0;
        for (DistributedVector::Iterator node_index = voltage.Begin();
             node_index!= voltage.End();
             ++node_index)
        {
            ChastePoint<2> location = monodomain_problem.rGetMesh().GetNode(node_index.Global)->GetPoint();

            if (fabs(location[0]-0.05)<1e-6 && fabs(location[1]-0.0)<1e-6) // Node 5
            {
                TS_ASSERT_DELTA(voltage[node_index.Global], 34.6692, 1e-3);
                edges_checked++;
            }

            if (fabs(location[0]-0.0)<1e-6 && fabs(location[1]-0.05)<1e-6) // Node 55
            {
                TS_ASSERT_DELTA(voltage[node_index.Global], 34.6692, 1e-3);
                edges_checked++;
            }

            if (fabs(location[0]-0.1)<1e-6 && fabs(location[1]-0.05)<1e-6) // Node 65
            {
                TS_ASSERT_DELTA(voltage[node_index.Global], 34.6692, 1e-3);
                edges_checked++;
            }

            if (fabs(location[0]-0.05)<1e-6 && fabs(location[1]-0.1)<1e-6) // Node 115
            {
                TS_ASSERT_DELTA(voltage[node_index.Global], 34.6692, 1e-3);
                edges_checked++;
            }
        }

        unsigned edges_checked_reduced;
        MPI_Allreduce(&edges_checked, &edges_checked_reduced, 1, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);
        TS_ASSERT(edges_checked_reduced==4);
    }

    // Same as TestMonodomainProblem1D, but uses NO matrix based assembly.
    void TestMonodomainWithNoMatrixBasedAssembly() throw(Exception)
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetPrintingTimeStep(1); //ms
        HeartConfig::Instance()->SetSimulationDuration(2); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("MonoProblem1dMatrixBasedRhs");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_1d");

        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1> cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );

        // switch off matrix based assembly
        monodomain_problem.UseMatrixBasedRhsAssembly(false);

        monodomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        monodomain_problem.Solve();

        // test whether voltages and gating variables are in correct ranges
        CheckMonoLr91Vars<1>(monodomain_problem);

        // check some voltages
        ReplicatableVector voltage_replicated(monodomain_problem.GetSolution());
        double atol=5e-3;

        TS_ASSERT_DELTA(voltage_replicated[1], 20.7710232, atol);
        TS_ASSERT_DELTA(voltage_replicated[3], 21.5319692, atol);
        TS_ASSERT_DELTA(voltage_replicated[5], 22.9280817, atol);
        TS_ASSERT_DELTA(voltage_replicated[7], 24.0611303, atol);
        TS_ASSERT_DELTA(voltage_replicated[9], -0.770330519, atol);
        TS_ASSERT_DELTA(voltage_replicated[10], -19.2234919, atol);
    }


    ///////////////////////////////////////////////////////////////////
    // Solve a simple simulation and check the output was only
    // printed out at the correct times
    ///////////////////////////////////////////////////////////////////
    void TestMonodomainProblemPrintsOnlyAtRequestedTimes() throw(Exception)
    {
        HeartConfig::Instance()->SetPrintingTimeStep(0.1);
        HeartConfig::Instance()->SetSimulationDuration(0.3); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("MonoProblem1d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("mono_testPrintTimes");

        // run testing PrintingTimeSteps
        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1> cell_factory;
        MonodomainProblem<1>* p_monodomain_problem = new MonodomainProblem<1>( &cell_factory );

        p_monodomain_problem->Initialise();
        p_monodomain_problem->SetWriteInfo();

        p_monodomain_problem->Solve();


        // read data entries for the time file and check correct
        //Hdf5DataReader data_reader1("MonoProblem1d", "mono_testPrintTimes");
        Hdf5DataReader data_reader1= p_monodomain_problem->GetDataReader();
        delete p_monodomain_problem;

        std::vector<double> times = data_reader1.GetUnlimitedDimensionValues();

        TS_ASSERT_EQUALS( times.size(), 4u);
        TS_ASSERT_DELTA( times[0], 0.00, 1e-12);
        TS_ASSERT_DELTA( times[1], 0.10, 1e-12);
        TS_ASSERT_DELTA( times[2], 0.20, 1e-12);
        TS_ASSERT_DELTA( times[3], 0.30, 1e-12);


    }

    void TestMonodomainWithMeshInMemoryToMeshalyzer() throw(Exception)
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005, 0.0005));
        HeartConfig::Instance()->SetSimulationDuration(0.1);  //ms
        HeartConfig::Instance()->SetOutputDirectory("Monodomain2d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("monodomain2d");
        
        //set the postprocessing information we want
        std::vector<unsigned> origin_nodes;
        origin_nodes.push_back(0u);
        HeartConfig::Instance()->SetConductionVelocityMaps(origin_nodes);

        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 2> cell_factory;

        ///////////////////////////////////////////////////////////////////
        // monodomain
        ///////////////////////////////////////////////////////////////////
        MonodomainProblem<2> monodomain_problem( &cell_factory );

        monodomain_problem.ConvertOutputToMeshalyzerFormat(true); // for coverage
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(10, 10, true);
        mesh.Scale(0.01,0.01); //To get 1mm x 1mm
        monodomain_problem.SetMesh(&mesh);
        monodomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        //Clean previous output
        OutputFileHandler handler("Monodomain2d/output", true);
        PetscTools::Barrier();

        // Checking that the following files don't exist in the output directory before calling Solve()
        unsigned num_files=6;
        std::string test_file_names[6]={"monodomain2d_mesh.pts", "monodomain2d_mesh.tri", "monodomain2d_V.dat",
              "ChasteParameters.xml", "ChasteDefaults.xml", "ConductionVelocityFromNode0.dat"};
        for (unsigned i=0; i<num_files; i++)
        {
            std::string compare_command = "cmp -s ";
            compare_command += handler.GetOutputDirectoryFullPath("Monodomain2d/output")+"/"+test_file_names[i];
            compare_command += " ";
            compare_command += "heart/test/data/Monodomain2d/";
            compare_command += test_file_names[i];
            TS_ASSERT_EQUALS(system(compare_command.c_str()), 512);//Not there
        }
        
        // now solve
        monodomain_problem.Solve();

        PetscTools::Barrier();
        // Compare output files
        for (unsigned i=0; i<num_files; i++)
        {
            if(test_file_names[i] == "monodomain2d_V.dat")
            {
                /*
                 * Since we started using bjacobi as the default preconditioner, parallel and sequential tests
                 * may return different results (always accurate to the tolerance requested). "diff" is unable
                 * to take this consideration into account.
                 *
                 * We will test that the file exists though.
                 */
                std::ifstream vm_file;
                std::string command = handler.GetOutputDirectoryFullPath("Monodomain2d/output")+"/"+test_file_names[i];
                vm_file.open(command.c_str());
                TS_ASSERT(vm_file.is_open());
                vm_file.close();
            }
            else
            {
                std::string compare_command = "diff --ignore-matching-lines=\"<ChasteParameters\" ";
                compare_command += handler.GetOutputDirectoryFullPath("Monodomain2d/output")+"/"+test_file_names[i];
                compare_command += " ";
                compare_command += "heart/test/data/Monodomain2d/";
                compare_command += test_file_names[i];
                TS_ASSERT_EQUALS(system(compare_command.c_str()), 0);
            }
        }
    }

    void TestMonodomainProblemCreates1DGeometry()
    {
        HeartConfig::Instance()->SetSpaceDimension(1);
        HeartConfig::Instance()->SetFibreLength(1, 0.01);
        HeartConfig::Instance()->SetSimulationDuration(5.0);  //ms
        HeartConfig::Instance()->SetOutputDirectory("MonodomainCreatesGeometry");
        HeartConfig::Instance()->SetOutputFilenamePrefix("monodomain1d");

        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1> cell_factory(-600 * 5000);

        MonodomainProblem<1> monodomain_problem( &cell_factory );

        monodomain_problem.ConvertOutputToMeshalyzerFormat(true);
        monodomain_problem.Initialise();

        monodomain_problem.Solve();

        // check some voltages
        ReplicatableVector voltage_replicated(monodomain_problem.GetSolution());
        double atol=5e-3;

        TS_ASSERT_DELTA(voltage_replicated[1], 13.9682, atol);
        TS_ASSERT_DELTA(voltage_replicated[3], 13.9149, atol);
        TS_ASSERT_DELTA(voltage_replicated[5], 13.8216, atol);
        TS_ASSERT_DELTA(voltage_replicated[7], 13.7275, atol);
    }

    void TestMonodomainProblemCreates2DGeometry()
    {
        HeartConfig::Instance()->SetSpaceDimension(2);
        HeartConfig::Instance()->SetSheetDimensions(0.3, 0.3, 0.01);
        HeartConfig::Instance()->SetSimulationDuration(4.0);  //ms
        HeartConfig::Instance()->SetOutputDirectory("MonodomainCreatesGeometry");
        HeartConfig::Instance()->SetOutputFilenamePrefix("monodomain2d");

        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 2> cell_factory(-600 * 5000);

        MonodomainProblem<2> monodomain_problem( &cell_factory );

        monodomain_problem.ConvertOutputToMeshalyzerFormat(true);
        monodomain_problem.Initialise();

        monodomain_problem.Solve();

        // check some voltages
        ReplicatableVector voltage_replicated(monodomain_problem.GetSolution());
        double atol=5e-3;

        TS_ASSERT_DELTA(voltage_replicated[1], 17.5728, atol);
        TS_ASSERT_DELTA(voltage_replicated[3], 17.4562, atol);
        TS_ASSERT_DELTA(voltage_replicated[5], 17.2559, atol);
        TS_ASSERT_DELTA(voltage_replicated[7], 17.0440, atol);
    }

    void TestMonodomainProblemCreates3DGeometry()
    {
        HeartConfig::Instance()->SetSpaceDimension(3);
        HeartConfig::Instance()->SetSlabDimensions(0.1, 0.1, 0.1, 0.01);
        HeartConfig::Instance()->SetSimulationDuration(2.0);  //ms
        HeartConfig::Instance()->SetOutputDirectory("MonodomainCreatesGeometry");
        HeartConfig::Instance()->SetOutputFilenamePrefix("monodomain3d");

        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 3> cell_factory(-600 * 5000);

        MonodomainProblem<3> monodomain_problem( &cell_factory );

        monodomain_problem.ConvertOutputToMeshalyzerFormat(true);
        monodomain_problem.Initialise();

        monodomain_problem.Solve();

        // check some voltages
        ReplicatableVector voltage_replicated(monodomain_problem.GetSolution());
        double atol=5e-3;

        TS_ASSERT_DELTA(voltage_replicated[1], 31.8629, atol);
        TS_ASSERT_DELTA(voltage_replicated[3], 32.0507, atol);
        TS_ASSERT_DELTA(voltage_replicated[5], 32.6340, atol);
        TS_ASSERT_DELTA(voltage_replicated[7], 33.6119, atol);
    }

    // Test the functionality for outputing the values of requested cell state variables
    void TestMonodomainProblemPrintsMultipleVariables() throw (Exception)
    {
        // Set configuration file 
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/MultipleVariablesMonodomain.xml");
   
        // Set up problem
        PlaneStimulusCellFactory<FaberRudy2000Version3, 1> cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );

        // Solve
        monodomain_problem.Initialise();
        monodomain_problem.Solve();
        monodomain_problem.ConvertOutputToMeshalyzerFormat();

        // Get a reference to a reader object for the simulation results
        Hdf5DataReader data_reader1=monodomain_problem.GetDataReader();
        std::vector<double> times = data_reader1.GetUnlimitedDimensionValues();

        // Check there is information about 101 timesteps (0, 0.01, 0.02, ...) 
        TS_ASSERT_EQUALS( times.size(), 11u);
        TS_ASSERT_DELTA( times[0], 0.0, 1e-12);
        TS_ASSERT_DELTA( times[1], 0.01, 1e-12);
        TS_ASSERT_DELTA( times[2], 0.02, 1e-12);
        TS_ASSERT_DELTA( times[3], 0.03, 1e-12);

        // There should be 101 values per variable and node.
        std::vector<double> node_5_v = data_reader1.GetVariableOverTime("V", 5);
        TS_ASSERT_EQUALS( node_5_v.size(), 11u);

        std::vector<double> node_5_cai = data_reader1.GetVariableOverTime("CaI", 5);
        TS_ASSERT_EQUALS( node_5_cai.size(), 11U);

        std::vector<double> node_5_nai = data_reader1.GetVariableOverTime("Nai", 5);
        TS_ASSERT_EQUALS( node_5_nai.size(), 11U);        

        std::vector<double> node_5_ki = data_reader1.GetVariableOverTime("Ki", 5);
        TS_ASSERT_EQUALS( node_5_ki.size(), 11U);
    }
    
    void TestMonodomainProblemExceptions() throw (Exception)
    {
        HeartConfig::Instance()->SetSimulationDuration(1.0); //ms

        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1> cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );

        // Throws because we've not called initialise
        TS_ASSERT_THROWS_THIS(monodomain_problem.Solve(),"Pde is null, Initialise() probably hasn\'t been called");

        // Throws because mesh filename is unset
        TS_ASSERT_THROWS_THIS(monodomain_problem.Initialise(),
                "No mesh given: define it in XML parameters file or call SetMesh()\n"
                "No Mesh provided (neither default nor user defined)");

        // Throws because initialise hasn't been called
        TS_ASSERT_THROWS_THIS(monodomain_problem.Solve(),"Pde is null, Initialise() probably hasn\'t been called");

        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("");
        HeartConfig::Instance()->SetOutputFilenamePrefix("");

        monodomain_problem.Initialise();

        //Throws because the HDF5 slab isn't on the disk
        TS_ASSERT_THROWS_THIS(monodomain_problem.GetDataReader(),"Data reader invalid as data writer cannot be initialised");

        // throw because end time is negative
        HeartConfig::Instance()->SetSimulationDuration(-1.0); //ms
        TS_ASSERT_THROWS_THIS(monodomain_problem.Solve(),"End time should be greater than 0");
        HeartConfig::Instance()->SetSimulationDuration( 1.0); //ms

        // throws because output dir and filename are both ""
        TS_ASSERT_THROWS_THIS(monodomain_problem.Solve(),
                "Either explicitly specify not to print output (call PrintOutput(false)) or "
                "specify the output directory and filename prefix");
    }
    
    /**
     * Not a very thorough test yet - just checks we can load a problem, simulate it, and
     * get expected results.
     */
    void TestArchiving() throw(Exception)
    {
        // Based on TestMonodomainProblem1D()
        OutputFileHandler handler("monadomain_problem_archive", false);
        handler.SetArchiveDirectory();
        std::string archive_filename = ArchiveLocationInfo::GetProcessUniqueFilePath("monodomain_problem.arch");

        // Save
        {
            HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
            HeartConfig::Instance()->SetSimulationDuration(2.0); //ms
            HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
            HeartConfig::Instance()->SetOutputDirectory("MonoProblemArchive");
            HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_1d");
    
            PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1> cell_factory;
            MonodomainProblem<1> monodomain_problem( &cell_factory );
    
            monodomain_problem.Initialise();
    
            HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
            HeartConfig::Instance()->SetCapacitance(1.0);
            
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            AbstractCardiacProblem<1,1,1>* const p_monodomain_problem = &monodomain_problem;
            output_arch & p_monodomain_problem;
        }
        
        // Load
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractCardiacProblem<1,1,1> *p_monodomain_problem;
            input_arch >> p_monodomain_problem;

            p_monodomain_problem->Solve();

            // test whether voltages and gating variables are in correct ranges
            CheckMonoLr91Vars<1>(static_cast<MonodomainProblem<1,1>&>(*p_monodomain_problem));
    
            // check some voltages
            ReplicatableVector voltage_replicated(p_monodomain_problem->GetSolution());
            double atol=5e-3;
    
            TS_ASSERT_DELTA(voltage_replicated[1], 20.7710232, atol);
            TS_ASSERT_DELTA(voltage_replicated[3], 21.5319692, atol);
            TS_ASSERT_DELTA(voltage_replicated[5], 22.9280817, atol);
            TS_ASSERT_DELTA(voltage_replicated[7], 24.0611303, atol);
            TS_ASSERT_DELTA(voltage_replicated[9], -0.770330519, atol);
            TS_ASSERT_DELTA(voltage_replicated[10], -19.2234919, atol);

            // check a progress report exists
            TS_ASSERT_EQUALS(system(("ls " + OutputFileHandler::GetChasteTestOutputDirectory() + "MonoProblemArchive/").c_str()), 0);
        }
    }
};

#endif //_TESTMONODOMAINPROBLEM_HPP_
