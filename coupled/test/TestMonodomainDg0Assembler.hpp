#ifndef _TESTMONODOMAINDG0ASSEMBLER_HPP_
#define _TESTMONODOMAINDG0ASSEMBLER_HPP_

// Element.hpp includes the Boost ublas objects - these need to
// be included early...  We think.  We're not that sure.
#include "Element.hpp"

#include <cxxtest/TestSuite.h>
#include <petsc.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <sys/stat.h>
#include <sys/types.h>

#include "SimpleLinearSolver.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "Node.hpp"
#include "Element.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "MonodomainDg0Assembler.hpp"
#include "TrianglesMeshReader.hpp"
#include "ColumnDataWriter.hpp"
#include "ColumnDataReader.hpp"
#include "PropagationPropertiesCalculator.hpp"

#include "MonodomainPde.hpp"
#include "MockEulerIvpOdeSolver.hpp"

#include "PetscSetupAndFinalize.hpp"
#include "MonodomainProblem.hpp"
#include "AbstractLinearParabolicPde.hpp"
#include "AbstractMonodomainProblemStimulus.hpp"

#include<time.h>
class PointStimulus1D: public AbstractMonodomainProblemStimulus<1>
{
public:
    virtual void Apply(MonodomainPde<1> *pPde, 
                       ConformingTetrahedralMesh<1,1> *pMesh)
    {
        static InitialStimulus stimulus(-600.0, 0.5);

        for (int i = 0; i < pMesh->GetNumNodes(); i++)
        {
            if (pMesh->GetNodeAt(i)->GetPoint()[0] == 0.0)
            {
                pPde->SetStimulusFunctionAtNode(i, &stimulus);
            }
        }
    }
    virtual ~PointStimulus1D()
    {
    }   
};

class EdgeStimulus2D: public AbstractMonodomainProblemStimulus<2>
{
public:
    virtual void Apply(MonodomainPde<2> *pPde,
                       ConformingTetrahedralMesh<2,2> *pMesh)
    {
        static InitialStimulus stimulus(-600.0, 0.5);
        
        for (int i = 0; i < pMesh->GetNumNodes(); i++)
        {
            if (pMesh->GetNodeAt(i)->GetPoint()[0] == 0.0)
            {
                pPde->SetStimulusFunctionAtNode(i, &stimulus);
            }
        }
    }
    virtual ~EdgeStimulus2D()
    {
    }
};

class PointStimulus2D: public AbstractMonodomainProblemStimulus<2>
{
private:
    int mNode;
    
public:
    PointStimulus2D(const int node = 0)
    {
        mNode = node;
    }
    
    virtual void Apply(MonodomainPde<2> *pPde,
                       ConformingTetrahedralMesh<2,2> *)
    {
        static InitialStimulus stimulus(-6000.0, 0.5);

        pPde->SetStimulusFunctionAtNode(mNode, &stimulus);
    }
    virtual ~PointStimulus2D()
    {
    }
};

class TestMonodomainDg0Assembler : public CxxTest::TestSuite 
{   
private:
	/**
	 * Refactor code to set up a PETSc vector holding the initial condition.
	 */
	Vec CreateInitialConditionVec(int size)
	{
    	Vec initial_condition;
    	VecCreate(PETSC_COMM_WORLD, &initial_condition);
    	VecSetSizes(initial_condition, PETSC_DECIDE, size);
    	VecSetFromOptions(initial_condition);
    	return initial_condition;
	}
    
public:

    // Solve on a 1D string of cells, 1mm long with a space step of 0.1mm.
    void TestMonodomainDg01D()
    {
        PointStimulus1D point_stimulus_1D;
        MonodomainProblem<1> monodomainProblem("mesh/test/data/1D_0_to_1mm_10_elements",
                                               2, // ms
                                               "testoutput/MonoDg01d",
                                               "NewMonodomainLR91_1d",
                                               &point_stimulus_1D);

        monodomainProblem.Solve();

        double* voltage_array;

        // test whether voltages and gating variables are in correct ranges

        VecGetArray(monodomainProblem.mCurrentVoltage, &voltage_array);

        for(int global_index=monodomainProblem.mLo; global_index<monodomainProblem.mHi; global_index++)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77
            double Ena   =  54.4;   // mV
            double Ek    = -77.0;   // mV

            TS_ASSERT_LESS_THAN_EQUALS(   voltage_array[global_index-monodomainProblem.mLo] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(  -voltage_array[global_index-monodomainProblem.mLo] + (Ek-30), 0);

            std::vector<double> odeVars = monodomainProblem.mMonodomainPde->GetOdeVarsAtNode(global_index);
            for(int j=0; j<8; j++)
            {
                // if not voltage or calcium ion conc, test whether between 0 and 1
                if((j!=4) && (j!=3))
                {
                    TS_ASSERT_LESS_THAN_EQUALS(  odeVars[j], 1.0);
                    TS_ASSERT_LESS_THAN_EQUALS( -odeVars[j], 0.0);
                }
            }

            if (global_index==1)
            {
                TS_ASSERT_DELTA(voltage_array[global_index-monodomainProblem.mLo], 20.0416, 0.001);
            }
            if (global_index==3)
            {
                TS_ASSERT_DELTA(voltage_array[global_index-monodomainProblem.mLo], 20.7888, 0.001);
            }
            if (global_index==5)
            {
                TS_ASSERT_DELTA(voltage_array[global_index-monodomainProblem.mLo], 22.0786, 0.001);
            }
            if (global_index==7)
            {
                TS_ASSERT_DELTA(voltage_array[global_index-monodomainProblem.mLo], 22.2901, 0.001);
            }
            if (global_index==9)
            {
                TS_ASSERT_DELTA(voltage_array[global_index-monodomainProblem.mLo], -15.7600, 0.001);
            }
            if (global_index==10) // RHS
            {
                TS_ASSERT_DELTA(voltage_array[global_index-monodomainProblem.mLo], -35.1363, 0.001);
            }
        }

        VecRestoreArray(monodomainProblem.mCurrentVoltage, &voltage_array);      
        VecAssemblyBegin(monodomainProblem.mCurrentVoltage);
        VecAssemblyEnd(monodomainProblem.mCurrentVoltage);
        VecDestroy(monodomainProblem.mCurrentVoltage);

    }
    
    // Solve on a 2D 1mm by 1mm mesh (space step = 0.1mm), stimulating the left
    // edge.
    // Should behave like the 1D case, extrapolated.
    void TestMonodomainDg02DWithEdgeStimulus( void )
    {   
        EdgeStimulus2D edge_stimulus_2D;
        
        // using the criss-cross mesh so wave propagates properly
        MonodomainProblem<2> monodomainProblem("mesh/test/data/2D_0_to_1mm_400_elements",
                                               2,   // ms
                                               "testoutput/MonoDg02dWithEdgeStimulus",
                                               "NewMonodomainLR91_2dWithEdgeStimulus",
                                               &edge_stimulus_2D);

        monodomainProblem.Solve();
        
        double* voltage_array;
        VecGetArray(monodomainProblem.mCurrentVoltage, &voltage_array); 
    
        // test whether voltages and gating variables are in correct ranges
        for(int global_index=monodomainProblem.mLo; global_index<monodomainProblem.mHi; global_index++)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77
            double Ena   =  54.4;
            double Ek    = -77.0;
            
            TS_ASSERT_LESS_THAN_EQUALS(   voltage_array[global_index-monodomainProblem.mLo] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(  -voltage_array[global_index-monodomainProblem.mLo] + (Ek-30), 0);
                
            std::vector<double> odeVars = monodomainProblem.mMonodomainPde->GetOdeVarsAtNode(global_index);           
            for(int j=0; j<8; j++)
            {
                // if not voltage or calcium ion conc, test whether between 0 and 1 
                if((j!=4) && (j!=3))
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
            for (int i = 0; i < monodomainProblem.mMesh.GetNumNodes(); i++)
            {
                if (monodomainProblem.mMesh.GetNodeAt(i)->GetPoint()[0] == 0.1)
                {
                    // x = 0 is where the stimulus has been applied
                    // x = 0.1cm is the other end of the mesh and where we want to 
                    //       to test the value of the nodes
                    
                    if (need_initialisation)
                    {
                        voltage = voltage_array[i];
                        need_initialisation = false;
                    }
                    else
                    {
                        // Tests the final voltages for all the RHS edge nodes
                        // are close to each other.
                        // This works as we are using the 'criss-cross' mesh,
                        // the voltages would vary more with a mesh with all the
                        // triangles aligned in the same direction.
                        TS_ASSERT_DELTA(voltage_array[i], voltage, 0.01);

                       // std::cout << "y=" << monodomainProblem.mMesh.GetNodeAt(i)->GetPoint()[1] << std::endl;
                    }
                    
                    
                    // Check against 1d case - THIS TEST HAS BEEN REMOVED AS THE MESH
                    // IS FINER THAN THE 1D MESH SO WE DONT EXPECT THE RESULTS TO BE THE SAME
                    // TS_ASSERT_DELTA(voltage_array[i], -35.1363, 35*0.1);
                    
                    // test the RHS edge voltages
                    // hardcoded result that looks accurate - this is a test to see
                    // that nothing has changeed
                    // assumes endtime = 2ms
                    TS_ASSERT_DELTA(voltage_array[i], -64.8086, 0.01);
                }
            }
        }
        
        VecRestoreArray(monodomainProblem.mCurrentVoltage, &voltage_array);      
        VecAssemblyBegin(monodomainProblem.mCurrentVoltage);
        VecAssemblyEnd(monodomainProblem.mCurrentVoltage);
        VecDestroy(monodomainProblem.mCurrentVoltage);
    }   

    // Solve on a 2D 1mm by 1mm mesh (space step = 0.05mm), stimulating the left
    // edge.
    // Should behave like the 1D case, extrapolated.
    void xTestMonodomainDg02DWithEdgeStimulusOnFinerMesh( void )
    {   
        EdgeStimulus2D edge_stimulus_2D;
        
        MonodomainProblem<2> monodomainProblem("mesh/test/data/2D_0_to_1mm_800_elements",
                                               2,   // ms
                                               "testoutput/MonoDg02dWithEdgeStimulus",
                                               "NewMonodomainLR91_2dWithEdgeStimulus",
                                               &edge_stimulus_2D);
        monodomainProblem.SetTimeSteps( 0.005, 0.005);
        monodomainProblem.Solve();
        
        double* voltage_array;
        VecGetArray(monodomainProblem.mCurrentVoltage, &voltage_array); 
    
        // test whether voltages and gating variables are in correct ranges
        for(int global_index=monodomainProblem.mLo; global_index<monodomainProblem.mHi; global_index++)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77
            double Ena   =  54.4;
            double Ek    = -77.0;
            
            TS_ASSERT_LESS_THAN_EQUALS(   voltage_array[global_index-monodomainProblem.mLo] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(  -voltage_array[global_index-monodomainProblem.mLo] + (Ek-30), 0);
                
            std::vector<double> odeVars = monodomainProblem.mMonodomainPde->GetOdeVarsAtNode(global_index);           
            for(int j=0; j<8; j++)
            {
                // if not voltage or calcium ion conc, test whether between 0 and 1 
                if((j!=4) && (j!=3))
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
            for (int i = 0; i < monodomainProblem.mMesh.GetNumNodes(); i++)
            {
                if (monodomainProblem.mMesh.GetNodeAt(i)->GetPoint()[0] == 0.1)
                {
                    // x = 0 is where the stimulus has been applied
                    // x = 0.1cm is the other end of the mesh and where we want to 
                    //       to test the value of the nodes
                    
                    if (need_initialisation)
                    {
                        voltage = voltage_array[i];
                        need_initialisation = false;
                    }
                    else
                    {
                        // Note that the RHS will be sampled during the upstroke,
                        // thus varies by about 3mV in 1 timestep.  So this isn't
                        // quite as bad as it looks.
                        // A finer mesh should give better results...
                        TS_ASSERT_DELTA(voltage_array[i], voltage, 4.0);
                       // std::cout << "y=" << monodomainProblem.mMesh.GetNodeAt(i)->GetPoint()[1] << std::endl;
                    }
                    
                    // Check against 1d case
                    TS_ASSERT_DELTA(voltage_array[i], -35.1363, 3.0);
                }
            }
        }
        
        VecRestoreArray(monodomainProblem.mCurrentVoltage, &voltage_array);      
        VecAssemblyBegin(monodomainProblem.mCurrentVoltage);
        VecAssemblyEnd(monodomainProblem.mCurrentVoltage);
        VecDestroy(monodomainProblem.mCurrentVoltage);
    }

    // Solve on a 2D 1mm by 1mm mesh (space step = 0.1mm), stimulating in the
    // very centre of the mesh.
    void TestMonodomainDg02DWithPointStimulusInTheVeryCentreOfTheMesh( void )
    {   
        // To time the solve
        time_t start,end;
        double dif;
        time (&start);
        
        
        PointStimulus2D point_stimulus_2D(60); // Central node
        
        MonodomainProblem<2> monodomainProblem("mesh/test/data/2D_0_to_1mm_400_elements",
                                               1.3,   // ms - needs to be 1.3 to pass test
                                               "testoutput/MonoDg02dWithPointStimulus",
                                               "NewMonodomainLR91_2dWithPointStimulus",
                                               &point_stimulus_2D);


        monodomainProblem.Solve();
        
        // To time the solve
        time (&end);
        dif = difftime (end,start);
        //printf ("\nSolve took %.2lf seconds. \n", dif );
        
        double* voltage_array;
    
        // test whether voltages and gating variables are in correct ranges

        VecGetArray(monodomainProblem.mCurrentVoltage, &voltage_array); 
        
        for(int global_index=monodomainProblem.mLo; global_index<monodomainProblem.mHi; global_index++)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77
            double Ena   =  54.4;
            double Ek    = -77.0;
            
            TS_ASSERT_LESS_THAN_EQUALS(   voltage_array[global_index-monodomainProblem.mLo] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(  -voltage_array[global_index-monodomainProblem.mLo] + (Ek-30), 0);
                
            std::vector<double> odeVars = monodomainProblem.mMonodomainPde->GetOdeVarsAtNode(global_index);           
            for(int j=0; j<8; j++)
            {
                // if not voltage or calcium ion conc, test whether between 0 and 1 
                if((j!=4) && (j!=3))
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
            TS_ASSERT_DELTA(voltage_array[0], voltage_array[10],  0.1);
            TS_ASSERT_DELTA(voltage_array[0], voltage_array[110], 0.1);
            TS_ASSERT_DELTA(voltage_array[0], voltage_array[120], 0.1);
            
            // centres of edges
            TS_ASSERT_DELTA(voltage_array[5], voltage_array[55],  0.1);
            TS_ASSERT_DELTA(voltage_array[5], voltage_array[65],  0.1);
            TS_ASSERT_DELTA(voltage_array[5], voltage_array[115], 0.1);
            
            // hardcoded result to check nothing has changed
            // assumes endtime = 1.3
            TS_ASSERT_DELTA(voltage_array[0], -45.301, 0.1);
                        
        }
        
        VecRestoreArray(monodomainProblem.mCurrentVoltage, &voltage_array);      
        VecAssemblyBegin(monodomainProblem.mCurrentVoltage);
        VecAssemblyEnd(monodomainProblem.mCurrentVoltage);
        VecDestroy(monodomainProblem.mCurrentVoltage);
    }   
};
#endif //_TESTMONODOMAINDG0ASSEMBLER_HPP_
