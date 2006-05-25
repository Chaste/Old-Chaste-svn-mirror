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
    
    AbstractCardiacCell* CreateCardiacCellForNode(int node)
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
    
    AbstractCardiacCell* CreateCardiacCellForNode(int node)
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
    int mNodeNum;
public:
    PointStimulus2dCellFactory(int nodeNum) : AbstractCardiacCellFactory<2>(0.01)
    {
        mpStimulus = new InitialStimulus(-6000.0, 0.5);
        mNodeNum = nodeNum;
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(int node)
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
        monodomain_problem.SetEndTime(2);   // 2 ms
        monodomain_problem.SetOutputDirectory("/tmp/testoutput/MonoDg01d");
        monodomain_problem.SetOutputFilenamePrefix("NewMonodomainLR91_1d");

        monodomain_problem.Initialise();
        monodomain_problem.Solve();

        double* p_voltage_array;
        int lo, hi;
        monodomain_problem.GetVoltageArray(&p_voltage_array, lo, hi); 
    
        // test whether voltages and gating variables are in correct ranges
        for (int global_index=lo; global_index<hi; global_index++)
        {
            int local_index = global_index - lo;
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
                TS_ASSERT_DELTA(p_voltage_array[local_index], 19.4178, 0.001);
            }
            if (global_index==3)
            {
                TS_ASSERT_DELTA(p_voltage_array[local_index], 20.2563, 0.001);
            }
            if (global_index==5)
            {
                TS_ASSERT_DELTA(p_voltage_array[local_index], 21.7166, 0.001);
            }
            if (global_index==7)
            {
                TS_ASSERT_DELTA(p_voltage_array[local_index], 22.3258, 0.001);
            }
            if (global_index==9)
            {
                TS_ASSERT_DELTA(p_voltage_array[local_index], -12.2253, 0.001);
            }
            if (global_index==10) // RHS
            {
                TS_ASSERT_DELTA(p_voltage_array[local_index], -31.4838, 0.001);
            }
        }

        monodomain_problem.RestoreVoltageArray(&p_voltage_array);
    }


    
    // Solve on a 2D 1mm by 1mm mesh (space step = 0.1mm), stimulating the left
    // edge.
    // Should behave like the 1D case, extrapolated.
    void TestMonodomainDg02DWithEdgeStimulus( void )
    {   
        EdgeStimulusCellFactory cell_factory;
        
        // using the criss-cross mesh so wave propagates properly
        MonodomainProblem<2> monodomain_problem( &cell_factory );

        monodomain_problem.SetMeshFilename("mesh/test/data/2D_0_to_1mm_400_elements");
        monodomain_problem.SetEndTime(2);   // 2 ms
        monodomain_problem.SetOutputDirectory("/tmp/testoutput/MonoDg02dWithEdgeStimulus");
        monodomain_problem.SetOutputFilenamePrefix("NewMonodomainLR91_2dWithEdgeStimulus");

        monodomain_problem.Initialise();
        
        monodomain_problem.Solve();
        
        double* p_voltage_array;
        int lo, hi;
        monodomain_problem.GetVoltageArray(&p_voltage_array, lo, hi); 
    
        // test whether voltages and gating variables are in correct ranges
        for (int global_index=lo; global_index<hi; global_index++)
        {
            int local_index = global_index - lo;
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

        int num_procs;
        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);

        if (num_procs == 1)
        {
            /*
             * Test the top right node against the right one in the 1D case, 
             * comparing voltage, and then test all the nodes on the right hand 
             * side of the square against the top right one, comparing voltage.
             */
            bool need_initialisation = true;
            double voltage;

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
                        voltage = p_voltage_array[i];
                        need_initialisation = false;
                    }
                    else
                    {
                        // Tests the final voltages for all the RHS edge nodes
                        // are close to each other.
                        // This works as we are using the 'criss-cross' mesh,
                        // the voltages would vary more with a mesh with all the
                        // triangles aligned in the same direction.
                        TS_ASSERT_DELTA(p_voltage_array[i], voltage, 0.01);

                       // std::cout << "y=" << monodomain_problem.mMesh.GetNodeAt(i)->GetPoint()[1] << std::endl;
                    }
                    
                    
                    // Check against 1d case - THIS TEST HAS BEEN REMOVED AS THE MESH
                    // IS FINER THAN THE 1D MESH SO WE DONT EXPECT THE RESULTS TO BE THE SAME
                    // TS_ASSERT_DELTA(p_voltage_array[i], -35.1363, 35*0.1);
                    
                    // test the RHS edge voltages
                    // hardcoded result that looks accurate - this is a test to see
                    // that nothing has changeed
                    // assumes endtime = 2ms
                    TS_ASSERT_DELTA(p_voltage_array[i], -63.1852, 0.01);
                }
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
        
        PointStimulus2dCellFactory cell_factory(60); // Central node
        
        MonodomainProblem<2> monodomain_problem( &cell_factory );

        monodomain_problem.SetMeshFilename("mesh/test/data/2D_0_to_1mm_400_elements");
        monodomain_problem.SetEndTime(1.3);   // 1.3 ms - needs to be 1.3 ms to pass test
        monodomain_problem.SetOutputDirectory("/tmp/testoutput/MonoDg02dWithPointStimulus");
        monodomain_problem.SetOutputFilenamePrefix("NewMonodomainLR91_2dWithPointStimulus");

        monodomain_problem.Initialise();
        monodomain_problem.Solve();
        
        // To time the solve
        time (&end);
        dif = difftime (end,start);
        //printf ("\nSolve took %.2lf seconds. \n", dif );
        
        double* p_voltage_array;
        int lo, hi;
        monodomain_problem.GetVoltageArray(&p_voltage_array, lo, hi); 
    
        // test whether voltages and gating variables are in correct ranges
        for (int global_index=lo; global_index<hi; global_index++)
        {
            int local_index = global_index - lo;
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
        
        int num_procs;
        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);

        if (num_procs == 1)
        {
            /*
             * Test that corners are 'equal', and centres of sides.
             * Irregularities in which way the triangles are oriented make
             * this rather difficult, especially since the edges are sampled
             * during the upstroke.
             */
             
            // corners
            TS_ASSERT_DELTA(p_voltage_array[0], p_voltage_array[10],  0.1);
            TS_ASSERT_DELTA(p_voltage_array[0], p_voltage_array[110], 0.1);
            TS_ASSERT_DELTA(p_voltage_array[0], p_voltage_array[120], 0.1);
            
            // centres of edges
            TS_ASSERT_DELTA(p_voltage_array[5], p_voltage_array[55],  0.1);
            TS_ASSERT_DELTA(p_voltage_array[5], p_voltage_array[65],  0.1);
            TS_ASSERT_DELTA(p_voltage_array[5], p_voltage_array[115], 0.1);
            
            // hardcoded result to check nothing has changed
            // assumes endtime = 1.3
            TS_ASSERT_DELTA(p_voltage_array[0], -43.1522, 0.1);
                        
        }
        
      monodomain_problem.RestoreVoltageArray(&p_voltage_array);
     }   
};

#endif //_TESTMONODOMAINDG0ASSEMBLER_HPP_
