#ifndef _TESTMONODOMAINDG0ASSEMBLER_HPP_
#define _TESTMONODOMAINDG0ASSEMBLER_HPP_


#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include <petscvec.h>
#include <vector>
#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "ColumnDataReader.hpp"
#include "ReplicatableVector.hpp"
#include "CheckMonoLr91Vars.hpp"
#include "PlaneStimulusCellFactory.hpp"

class PointStimulus2dCellFactory : public AbstractCardiacCellFactory<2>
{
private:
    InitialStimulus *mpStimulus;
    unsigned mNodeNum;
public:
    PointStimulus2dCellFactory(int nodeNum) : AbstractCardiacCellFactory<2>(0.01)
    {
        mpStimulus = new InitialStimulus(-6000.0, 0.5);
        mNodeNum = nodeNum;
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        if (node == mNodeNum)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpStimulus);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpZeroStimulus);
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

    // Solve on a 1D string of cells, 1mm long with a space step of 0.1mm.
    void TestMonodomainDg01D()
    {
        PlaneStimulusCellFactory<1> cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );
        
        monodomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1mm_10_elements");
        monodomain_problem.SetEndTime(2);   // ms
        monodomain_problem.SetOutputDirectory("MonoDg01d");
        monodomain_problem.SetOutputFilenamePrefix("NewMonodomainLR91_1d");
        
        monodomain_problem.Initialise();
        
        monodomain_problem.GetMonodomainPde()->SetSurfaceAreaToVolumeRatio(1.0);
        monodomain_problem.GetMonodomainPde()->SetCapacitance(1.0);
        monodomain_problem.GetMonodomainPde()->SetIntracellularConductivityTensor(0.0005*identity_matrix<double>(1));
        
        monodomain_problem.Solve();
        
        // test whether voltages and gating variables are in correct ranges
        CheckMonoLr91Vars<1>(monodomain_problem);
        
        // check some voltages    
        ReplicatableVector voltage_replicated(monodomain_problem.GetVoltage());
        TS_ASSERT_DELTA(voltage_replicated[1], 20.7709, 0.001);
        TS_ASSERT_DELTA(voltage_replicated[3], 21.5321, 0.001);
        TS_ASSERT_DELTA(voltage_replicated[5], 22.9282, 0.001);
        TS_ASSERT_DELTA(voltage_replicated[7], 24.0612, 0.001);
        TS_ASSERT_DELTA(voltage_replicated[9], -0.7694, 0.001);
        TS_ASSERT_DELTA(voltage_replicated[10], -19.2224, 0.001);
       
    }
    
    // Solve on a 2D 1mm by 1mm mesh (space step = 0.1mm), stimulating the left
    // edge.
    // Should behave like the 1D case, extrapolated.
    // See also TestMonodomainSlab.hpp (nightly test) for the 3D version.
    void TestMonodomainDg02DWithEdgeStimulus( void )
    {
        static double test_tolerance=1e-10;
        PlaneStimulusCellFactory<2> cell_factory;
        
        // using the criss-cross mesh so wave propagates properly
        MonodomainProblem<2> monodomain_problem( &cell_factory );
        
        monodomain_problem.SetMeshFilename("mesh/test/data/2D_0_to_1mm_400_elements");
        monodomain_problem.SetEndTime(2);   // 2 ms
        monodomain_problem.SetOutputDirectory("MonoDg02dWithEdgeStimulus");
        monodomain_problem.SetOutputFilenamePrefix("NewMonodomainLR91_2dWithEdgeStimulus");
        
        monodomain_problem.Initialise();
        
        monodomain_problem.GetMonodomainPde()->SetSurfaceAreaToVolumeRatio(1.0);
        monodomain_problem.GetMonodomainPde()->SetCapacitance(1.0);
        monodomain_problem.GetMonodomainPde()->SetIntracellularConductivityTensor(0.0005*identity_matrix<double>(2));
        
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
        double probe_voltage;
        
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
                TS_ASSERT_DELTA(voltage_replicated[i], -59.6495, 1e-4);
            }
        }
        
    }
    
    
    // Solve on a 2D 1mm by 1mm mesh (space step = 0.1mm), stimulating in the
    // very centre of the mesh.
    void TestMonodomainDg02DWithPointStimulusInTheVeryCentreOfTheMesh( void )
    {
        // To time the solve
        time_t start,end;
        double dif;
        time (&start);
        static double test_tolerance=1e-10;
        
        PointStimulus2dCellFactory cell_factory(60); // Central node
        
        MonodomainProblem<2> monodomain_problem( &cell_factory );
        
        monodomain_problem.SetMeshFilename("mesh/test/data/2D_0_to_1mm_400_elements");
        monodomain_problem.SetEndTime(1.3);   // 1.3 ms - needs to be 1.3 ms to pass test
        monodomain_problem.SetOutputDirectory("MonoDg02dWithPointStimulus");
        monodomain_problem.SetOutputFilenamePrefix("NewMonodomainLR91_2dWithPointStimulus");
        
        monodomain_problem.Initialise();
        
        monodomain_problem.GetMonodomainPde()->SetSurfaceAreaToVolumeRatio(1.0);
        monodomain_problem.GetMonodomainPde()->SetCapacitance(1.0);
        monodomain_problem.GetMonodomainPde()->SetIntracellularConductivityTensor(0.0005*identity_matrix<double>(2));
        
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
        TS_ASSERT_DELTA(voltage_replicated[0], -34.3493, 1e-4);
        
//        monodomain_problem.RestoreVoltageArray(&p_voltage_array);
    }
    
    
    ///////////////////////////////////////////////////////////////////
    // Solve a simple simulation and check the output was only
    // printed out at the correct times
    ///////////////////////////////////////////////////////////////////
    void TestMonodomainProblemPrintsOnlyAtRequestedTimes()
    {
        // run testing PrintingTimeSteps
        PlaneStimulusCellFactory<1> cell_factory;
        MonodomainProblem<1>* p_monodomain_problem = new MonodomainProblem<1>( &cell_factory );
        
        p_monodomain_problem->SetMeshFilename("mesh/test/data/1D_0_to_1mm_10_elements");
        
        p_monodomain_problem->SetEndTime(0.30);          // ms
        p_monodomain_problem->SetPdeTimeStep(0.01);      // ms
        p_monodomain_problem->SetPrintingTimeStep(0.1);  // every 0.1ms
        
        p_monodomain_problem->SetOutputDirectory("MonoDg01d");
        p_monodomain_problem->SetOutputFilenamePrefix("mono_testPrintTimes");
        
        p_monodomain_problem->Initialise();
        p_monodomain_problem->Solve();
        
        delete p_monodomain_problem;
        
        // read data entries for the time file and check correct
        ColumnDataReader data_reader1("MonoDg01d", "mono_testPrintTimes");
        std::vector<double> times = data_reader1.GetUnlimitedDimensionValues();
        
        TS_ASSERT_EQUALS( times.size(), 4u);
        TS_ASSERT_DELTA( times[0], 0.00, 1e-12);
        TS_ASSERT_DELTA( times[1], 0.10, 1e-12);
        TS_ASSERT_DELTA( times[2], 0.20, 1e-12);
        TS_ASSERT_DELTA( times[3], 0.30, 1e-12);
        
        // run testing PrintEveryNthTimeStep
        p_monodomain_problem = new MonodomainProblem<1>( &cell_factory );
        
        p_monodomain_problem->SetMeshFilename("mesh/test/data/1D_0_to_1mm_10_elements");
        p_monodomain_problem->SetEndTime(0.50);   // ms
        p_monodomain_problem->SetOutputDirectory("MonoDg01d");
        p_monodomain_problem->SetOutputFilenamePrefix("mono_testPrintTimes");
        
        p_monodomain_problem->SetPdeTimeStep(0.01);
        p_monodomain_problem->PrintEveryNthTimeStep(17);  // every 17 timesteps
        
        p_monodomain_problem->SetWriteInfo(); // just to have SetWriteInfo() covered in the tests
        
        p_monodomain_problem->Initialise();
        p_monodomain_problem->Solve();
        
        // read data entries for the time file and check correct
        ColumnDataReader data_reader2("MonoDg01d", "mono_testPrintTimes");
        times = data_reader2.GetUnlimitedDimensionValues();
        
        TS_ASSERT_EQUALS( times.size(), 4u);
        TS_ASSERT_DELTA( times[0], 0.00,  1e-12);
        TS_ASSERT_DELTA( times[1], 0.17,  1e-12);
        TS_ASSERT_DELTA( times[2], 0.34,  1e-12);
        TS_ASSERT_DELTA( times[3], 0.50,  1e-12);
        
        delete p_monodomain_problem;
    }
    
    
    void TestMonodomainProblemExceptions() throw (Exception)
    {
        PlaneStimulusCellFactory<1> cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );

        // bad params
        TS_ASSERT_THROWS_ANYTHING(monodomain_problem.GetMonodomainPde()->SetSurfaceAreaToVolumeRatio(-1));
        
        TS_ASSERT_THROWS_ANYTHING(monodomain_problem.GetMonodomainPde()->SetCapacitance(-1));
        
        //Throws because we've not called initialise
        TS_ASSERT_THROWS_ANYTHING(monodomain_problem.Solve());
        
        // throws because argument is negative
        TS_ASSERT_THROWS_ANYTHING(monodomain_problem.SetPdeTimeStep(-1));
        
        // throws because argument is negative
        TS_ASSERT_THROWS_ANYTHING(monodomain_problem.SetPrintingTimeStep(-1));
        
        //Throws because mesh filename is unset
        TS_ASSERT_THROWS_ANYTHING(monodomain_problem.Initialise());
        monodomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1mm_10_elements");
        TS_ASSERT_THROWS_NOTHING(monodomain_problem.Initialise());
        
        //Throws because EndTime has not been set
        TS_ASSERT_THROWS_ANYTHING(monodomain_problem.Solve());
        monodomain_problem.SetEndTime(1);  // ms
    }
    
};

#endif //_TESTMONODOMAINDG0ASSEMBLER_HPP_
