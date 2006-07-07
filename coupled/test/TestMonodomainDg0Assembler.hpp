#ifndef _TESTMONODOMAINDG0ASSEMBLER_HPP_
#define _TESTMONODOMAINDG0ASSEMBLER_HPP_

// Element.hpp includes the Boost ublas objects - these need to
// be included early...  We think.  We're not that sure.
#include "Element.hpp"
#include <cxxtest/TestSuite.h>
#include <petscvec.h>
#include <vector>
//#include <iostream>
//#include "ConformingTetrahedralMesh.cpp"
#include "PetscSetupAndFinalize.hpp"
#include "MonodomainProblem.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "ColumnDataReader.hpp"
#include "ReplicatableVector.hpp" 

#include <time.h>

class PointStimulusCellFactory : public AbstractCardiacCellFactory<1>
{
private:
    // define a new stimulus
    InitialStimulus* mpStimulus;
    
public:
    PointStimulusCellFactory() : AbstractCardiacCellFactory<1>(0.01)
    {
        // set the new stimulus
        mpStimulus = new InitialStimulus(-600, 0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        if (mpMesh->GetNodeAt(node)->GetPoint()[0] == 0.0)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpStimulus);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpZeroStimulus);
        }
    }
    
    ~PointStimulusCellFactory(void)
    {
        delete mpStimulus;
    }
};


class EdgeStimulusCellFactory : public AbstractCardiacCellFactory<2>
{
private:
    InitialStimulus *mpStimulus;
public:
    EdgeStimulusCellFactory() : AbstractCardiacCellFactory<2>(0.01)
    {
        mpStimulus = new InitialStimulus(-600.0, 0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        if (mpMesh->GetNodeAt(node)->GetPoint()[0] == 0.0)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpStimulus);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpZeroStimulus);
        }
    }
    
    ~EdgeStimulusCellFactory(void)
    {
        delete mpStimulus;
    }
};

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
        PointStimulusCellFactory cell_factory;
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

        double* p_voltage_array;
        unsigned lo, hi;
        monodomain_problem.GetVoltageArray(&p_voltage_array, lo, hi); 
    
        // test whether voltages and gating variables are in correct ranges
        for (unsigned global_index=lo; global_index<hi; global_index++)
        {
            unsigned local_index = global_index - lo;
            // assuming LR model has Ena = 54.4 and Ek = -77
            double Ena   =  54.4;   // mV
            double Ek    = -77.0;   // mV

            TS_ASSERT_LESS_THAN_EQUALS( p_voltage_array[local_index] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(-p_voltage_array[local_index] + (Ek-30), 0);

            std::vector<double> odeVars = monodomain_problem.GetMonodomainPde()->GetCardiacCell(global_index)->rGetStateVariables();
            for (int j=0; j<8; j++)
            {
                // if not voltage or calcium ion conc, test whether between 0 and 1
                if ((j!=4) && (j!=3))
                {
                    TS_ASSERT_LESS_THAN_EQUALS(  odeVars[j], 1.0);
                    TS_ASSERT_LESS_THAN_EQUALS( -odeVars[j], 0.0);
                }
            }

            if (global_index==1)
            {
                TS_ASSERT_DELTA(p_voltage_array[local_index], 20.6690, 0.001);
            }
            if (global_index==3)
            {
                TS_ASSERT_DELTA(p_voltage_array[local_index], 21.4655, 0.001);
            }
            if (global_index==5)
            {
                TS_ASSERT_DELTA(p_voltage_array[local_index], 22.9016, 0.001);
            }
            if (global_index==7)
            {
                TS_ASSERT_DELTA(p_voltage_array[local_index], 24.0518, 0.001);
            }
            if (global_index==9)
            {
                TS_ASSERT_DELTA(p_voltage_array[local_index], -0.9282, 0.001);
            }
            if (global_index==10) // RHS
            {
                TS_ASSERT_DELTA(p_voltage_array[local_index], -19.4217, 0.001);
            }
        }
        monodomain_problem.RestoreVoltageArray(&p_voltage_array);
    }
    
    // Solve on a 2D 1mm by 1mm mesh (space step = 0.1mm), stimulating the left
    // edge.
    // Should behave like the 1D case, extrapolated.
    // See also TestMonodomainSlab.hpp (nightly test) for the 3D version.
    void TestMonodomainDg02DWithEdgeStimulus( void )
    {   
        static double test_tolerance=1e-10;
        EdgeStimulusCellFactory cell_factory;
        
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
        
        double* p_voltage_array;
        unsigned lo, hi;
        monodomain_problem.GetVoltageArray(&p_voltage_array, lo, hi); 
    
        // test whether voltages and gating variables are in correct ranges
        for (unsigned global_index=lo; global_index<hi; global_index++)
        {
            unsigned local_index = global_index - lo;
            // assuming LR model has Ena = 54.4 and Ek = -77
            double Ena   =  54.4;
            double Ek    = -77.0;
            
            TS_ASSERT_LESS_THAN_EQUALS( p_voltage_array[local_index] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(-p_voltage_array[local_index] + (Ek-30), 0);
                
            std::vector<double> odeVars = monodomain_problem.GetMonodomainPde()->GetCardiacCell(global_index)->rGetStateVariables();
            for (int j=0; j<8; j++)
            {
                // if not voltage or calcium ion conc, test whether between 0 and 1 
                if ((j!=4) && (j!=3))
                {
                    TS_ASSERT_LESS_THAN_EQUALS(  odeVars[j], 1.0);        
                    TS_ASSERT_LESS_THAN_EQUALS( -odeVars[j], 0.0);        
                }
            }
        }
        

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
        for (int i = 0; i < monodomain_problem.rGetMesh().GetNumNodes(); i++)
        {
            if (monodomain_problem.rGetMesh().GetNodeAt(i)->GetPoint()[0] == 0.1)
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
                TS_ASSERT_DELTA(voltage_replicated[i], -59.7978, 1e-4);
            }
        }
        
        monodomain_problem.RestoreVoltageArray(&p_voltage_array);
        
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
        
        double* p_voltage_array;
        unsigned lo, hi;
        monodomain_problem.GetVoltageArray(&p_voltage_array, lo, hi); 
    
        // test whether voltages and gating variables are in correct ranges
        for (unsigned global_index=lo; global_index<hi; global_index++)
        {
            unsigned local_index = global_index - lo;
            // assuming LR model has Ena = 54.4 and Ek = -77
            double Ena   =  54.4;
            double Ek    = -77.0;
            
            TS_ASSERT_LESS_THAN_EQUALS( p_voltage_array[local_index] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(-p_voltage_array[local_index] + (Ek-30), 0);
                
            std::vector<double> odeVars = monodomain_problem.GetMonodomainPde()->GetCardiacCell(global_index)->rGetStateVariables();
            for (int j=0; j<8; j++)
            {
                // if not voltage or calcium ion conc, test whether between 0 and 1 
                if ((j!=4) && (j!=3))
                {
                    TS_ASSERT_LESS_THAN_EQUALS(  odeVars[j], 1.0);        
                    TS_ASSERT_LESS_THAN_EQUALS( -odeVars[j], 0.0);        
                }
            }
        }
        
        

        /*
         * Test that corners are 'equal', and centres of sides.
         * Irregularities in which way the triangles are oriented make
         * this rather difficult, especially since the edges are sampled
         * during the upstroke.
         */
        Vec voltage=monodomain_problem.GetVoltage(); 
        ReplicatableVector voltage_replicated; 
        voltage_replicated.ReplicatePetscVector(voltage); 
          
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
        TS_ASSERT_DELTA(voltage_replicated[0], -34.7497, 1e-4);
                        
        monodomain_problem.RestoreVoltageArray(&p_voltage_array);
    }
    
    
    ///////////////////////////////////////////////////////////////////
    // Solve a simple simulation and check the output was only
    // printed out at the correct times
    ///////////////////////////////////////////////////////////////////
    void TestMonodomainProblemPrintsOnlyAtRequestedTimes()
    {
        // run testing PrintingTimeSteps
        PointStimulusCellFactory cell_factory;
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
        
        TS_ASSERT_EQUALS( times.size(), 4);
        TS_ASSERT_DELTA( times[0], 0.00, 1e-12);
        TS_ASSERT_DELTA( times[1], 0.10, 1e-12);
        TS_ASSERT_DELTA( times[2], 0.20, 1e-12);
        TS_ASSERT_DELTA( times[3], 0.30, 1e-12);

        // run testing PrintEveryNthTimeStep
        p_monodomain_problem = new MonodomainProblem<1>( &cell_factory );

        p_monodomain_problem->SetMeshFilename("mesh/test/data/1D_0_to_1mm_10_elements");
        p_monodomain_problem->SetEndTime(0.51);   // ms
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
                  
        TS_ASSERT_EQUALS( times.size(), 4);
        TS_ASSERT_DELTA( times[0], 0.00,  1e-12);
        TS_ASSERT_DELTA( times[1], 0.17,  1e-12);
        TS_ASSERT_DELTA( times[2], 0.34,  1e-12);
        TS_ASSERT_DELTA( times[3], 0.51,  1e-12);
        
        delete p_monodomain_problem;
    }        
};

#endif //_TESTMONODOMAINDG0ASSEMBLER_HPP_
