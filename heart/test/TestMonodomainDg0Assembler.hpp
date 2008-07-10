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


#ifndef _TESTMONODOMAINDG0ASSEMBLER_HPP_
#define _TESTMONODOMAINDG0ASSEMBLER_HPP_


#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include <petscvec.h>
#include <vector>
#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "Hdf5DataReader.hpp"
#include "ReplicatableVector.hpp"
#include "CheckMonoLr91Vars.hpp"
#include "PlaneStimulusCellFactory.hpp"


class PointStimulus2dCellFactory : public AbstractCardiacCellFactory<2>
{
private:
    SimpleStimulus *mpStimulus;
    unsigned mNodeNum;
public:
    PointStimulus2dCellFactory(int nodeNum) : AbstractCardiacCellFactory<2>()
    {
        mpStimulus = new SimpleStimulus(-6000.0, 0.5);
        mNodeNum = nodeNum;
    }

    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        if (node == mNodeNum)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpStimulus);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpZeroStimulus);
        }
    }

    ~PointStimulus2dCellFactory(void)
    {
        delete mpStimulus;
    }
};



class TestMonodomainDg0Assembler : public CxxTest::TestSuite
{
public:
    void tearDown()
    {
        HeartConfig::Destroy();   
    }

    // Solve on a 1D string of cells, 1mm long with a space step of 0.1mm.
    void TestMonodomainDg01D()
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetSimulationDuration(2.0); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");        
        HeartConfig::Instance()->SetOutputDirectory("MonoDg01d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_1d");
        
        PlaneStimulusCellFactory<1> cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );

        monodomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        monodomain_problem.Solve();

        // test whether voltages and gating variables are in correct ranges
        CheckMonoLr91Vars<1>(monodomain_problem);

        // check some voltages
        ReplicatableVector voltage_replicated(monodomain_problem.GetVoltage());
        double atol=5e-3;

        TS_ASSERT_DELTA(voltage_replicated[1], 20.7710232, atol);
        TS_ASSERT_DELTA(voltage_replicated[3], 21.5319692, atol);
        TS_ASSERT_DELTA(voltage_replicated[5], 22.9280817, atol);
        TS_ASSERT_DELTA(voltage_replicated[7], 24.0611303, atol);
        TS_ASSERT_DELTA(voltage_replicated[9], -0.770330519, atol);
        TS_ASSERT_DELTA(voltage_replicated[10], -19.2234919, atol);

        // cover get pde
        monodomain_problem.GetPde();
    }

    void TestMonodomainDg01DWithRelativeTolerance()
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetSimulationDuration(2.0); //ms
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);
        HeartConfig::Instance()->SetUseRelativeTolerance(1e-9);
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("MonoDg01d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_1d");
        
        PlaneStimulusCellFactory<1> cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );

        monodomain_problem.Initialise();

        monodomain_problem.Solve();

        // test whether voltages and gating variables are in correct ranges
        CheckMonoLr91Vars<1>(monodomain_problem);

        // check some voltages
        ReplicatableVector voltage_replicated(monodomain_problem.GetVoltage());
        double atol=1e-6;

        TS_ASSERT_DELTA(voltage_replicated[1], 20.7710232, atol);
        TS_ASSERT_DELTA(voltage_replicated[3], 21.5319692, atol);
        TS_ASSERT_DELTA(voltage_replicated[5], 22.9280817, atol);
        TS_ASSERT_DELTA(voltage_replicated[7], 24.0611303, atol);
        TS_ASSERT_DELTA(voltage_replicated[9], -0.770330519, atol);
        TS_ASSERT_DELTA(voltage_replicated[10], -19.2234919, atol);

    }

    void TestMonodomainDg01DWithAbsoluteTolerance() throw (Exception)
    {
        double atol = 1e-1;
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetSimulationDuration(2); //ms
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);
        HeartConfig::Instance()->SetUseAbsoluteTolerance(atol/4.0);
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("MonoDg01d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_1d");
        
        PlaneStimulusCellFactory<1> cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );

        monodomain_problem.Initialise();

        monodomain_problem.Solve();

        // test whether voltages and gating variables are in correct ranges
        CheckMonoLr91Vars<1>(monodomain_problem);

        // check some voltages
        ReplicatableVector voltage_replicated(monodomain_problem.GetVoltage());

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
    void TestMonodomainDg02DWithEdgeStimulus( void )
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005, 0.0005));
        HeartConfig::Instance()->SetSimulationDuration(2); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/2D_0_to_1mm_400_elements");
        HeartConfig::Instance()->SetOutputDirectory("MonoDg02dWithEdgeStimulus");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_2dWithEdgeStimulus");
       
        static double test_tolerance=1e-10;
        PlaneStimulusCellFactory<2> cell_factory;

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
        Vec voltage=monodomain_problem.GetVoltage();
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
        for (unsigned i = 0; i < monodomain_problem.rGetMesh().GetNumNodes(); i++)
        {
            if (monodomain_problem.rGetMesh().GetNode(i)->GetPoint()[0] == 0.1)
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
                    TS_ASSERT_DELTA(voltage_replicated[i], probe_voltage, test_tolerance);
                }


                // Check against 1d case - THIS TEST HAS BEEN REMOVED AS THE MESH
                // IS FINER THAN THE 1D MESH SO WE DONT EXPECT THE RESULTS TO BE THE SAME
                // TS_ASSERT_DELTA(p_voltage_array[i], -35.1363, 35*0.1);

                // test the RHS edge voltages
                // hardcoded result that looks accurate - this is a test to see
                // that nothing has changeed
                // assumes endtime = 2ms
                TS_ASSERT_DELTA(voltage_replicated[i], -59.6495, 4e-4);
            }
        }

    }


    // Solve on a 2D 1mm by 1mm mesh (space step = 0.1mm), stimulating in the
    // very centre of the mesh.
    void TestMonodomainDg02DWithPointStimulusInTheVeryCentreOfTheMesh( void )
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005, 0.0005));
        HeartConfig::Instance()->SetSimulationDuration(1.3); //ms - needs to be 1.3 ms to pass test
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/2D_0_to_1mm_400_elements");
        HeartConfig::Instance()->SetOutputDirectory("MonoDg02dWithPointStimulus");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_2dWithPointStimulus");
        
        // To time the solve
        time_t start,end;
        double dif;
        time (&start);
        static double test_tolerance=1e-10;

        PointStimulus2dCellFactory cell_factory(60); // Central node

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
        ReplicatableVector voltage_replicated(monodomain_problem.GetVoltage());

        // corners
        TS_ASSERT_DELTA(voltage_replicated[0], voltage_replicated[10],  test_tolerance);
        TS_ASSERT_DELTA(voltage_replicated[0], voltage_replicated[110], test_tolerance);
        TS_ASSERT_DELTA(voltage_replicated[0], voltage_replicated[120], test_tolerance);

        // centres of edges
        TS_ASSERT_DELTA(voltage_replicated[5], voltage_replicated[55],  test_tolerance);
        TS_ASSERT_DELTA(voltage_replicated[5], voltage_replicated[65],  test_tolerance);
        TS_ASSERT_DELTA(voltage_replicated[5], voltage_replicated[115], test_tolerance);

        // hardcoded result to check nothing has changed
        // assumes endtime = 1.3
        TS_ASSERT_DELTA(voltage_replicated[0], -34.3493, 1e-3);

//        monodomain_problem.RestoreVoltageArray(&p_voltage_array);
    }


    ///////////////////////////////////////////////////////////////////
    // Solve a simple simulation and check the output was only
    // printed out at the correct times
    ///////////////////////////////////////////////////////////////////
    void TestMonodomainProblemPrintsOnlyAtRequestedTimes()
    {
        HeartConfig::Instance()->SetPrintingTimeStep(0.1);
        HeartConfig::Instance()->SetSimulationDuration(0.3); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("MonoDg01d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("mono_testPrintTimes");
        
        // run testing PrintingTimeSteps
        PlaneStimulusCellFactory<1> cell_factory;
        MonodomainProblem<1>* p_monodomain_problem = new MonodomainProblem<1>( &cell_factory );

//        // pass in mesh this time rather than a filename
//        TrianglesMeshReader<1,1> reader("mesh/test/data/1D_0_to_1mm_10_elements");
//        ConformingTetrahedralMesh<1,1> mesh;
//        mesh.ConstructFromMeshReader(reader);
//        p_monodomain_problem->SetMesh(&mesh);

        p_monodomain_problem->Initialise();
        p_monodomain_problem->SetWriteInfo();
        
        p_monodomain_problem->Solve();

        delete p_monodomain_problem;

        // read data entries for the time file and check correct
        Hdf5DataReader data_reader1("MonoDg01d", "mono_testPrintTimes");
        std::vector<double> times = data_reader1.GetUnlimitedDimensionValues();

        TS_ASSERT_EQUALS( times.size(), 4u);
        TS_ASSERT_DELTA( times[0], 0.00, 1e-12);
        TS_ASSERT_DELTA( times[1], 0.10, 1e-12);
        TS_ASSERT_DELTA( times[2], 0.20, 1e-12);
        TS_ASSERT_DELTA( times[3], 0.30, 1e-12);

    }


    void TestMonodomainProblemExceptions() throw (Exception)
    {
        HeartConfig::Instance()->SetSimulationDuration(1.0); //ms
        
        PlaneStimulusCellFactory<1> cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );

        // Throws because we've not called initialise
        TS_ASSERT_THROWS_ANYTHING(monodomain_problem.Solve());

        // Throws because mesh filename is unset
        TS_ASSERT_THROWS_ANYTHING(monodomain_problem.Initialise());

        // Throws because output directory has not been set
        TS_ASSERT_THROWS_ANYTHING(monodomain_problem.Solve());
    }
};

#endif //_TESTMONODOMAINDG0ASSEMBLER_HPP_
