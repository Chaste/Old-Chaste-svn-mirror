#ifndef _TESTMONODOMAINDG0ASSEMBLERLONG_HPP_
#define _TESTMONODOMAINDG0ASSEMBLERLONG_HPP_

// Element.hpp includes the Boost ublas objects - these need to
// be included early...  We think.  We're not that sure.
#include "Element.hpp"

#include <cxxtest/TestSuite.h>
#include "petscvec.h"
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
                       ConformingTetrahedralMesh<2,2> *pMesh)
    {
        static InitialStimulus stimulus(-6000.0, 0.5);

        pPde->SetStimulusFunctionAtNode(mNode, &stimulus);
    }
};

class TestMonodomainDg0AssemblerLong : public CxxTest::TestSuite 
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


    // Solve on a 2D 1mm by 1mm mesh (space step = 0.1mm), stimulating in the
    // very centre of the mesh.
    // We run for 500 ms and then check that all the voltages at the final time
    // have returned to the resting potential of -84.5
    // test should take about 30mins (or less)
    void TestMonodomainDg02DWithPointStimulusInTheVeryCentreOfTheMesh( void )
    {   
        // To time the solve
        time_t start,end;
        double dif;
        time (&start);
        
        
        PointStimulus2D point_stimulus_2D(60); // Central node

        MonodomainProblem<2> monodomainProblem;

        monodomainProblem.SetMeshFilename("mesh/test/data/2D_0_to_1mm_400_elements");
        monodomainProblem.SetEndTime(500);   // 500 ms
        monodomainProblem.SetOutputDirectory("testoutput/MonoDg02dWithPointStimulusLong");
        monodomainProblem.SetOutputFilenamePrefix("NewMonodomainLR91_2dWithPointStimulusLong");
        monodomainProblem.SetStimulus(&point_stimulus_2D);

        monodomainProblem.Solve();
        
        // To time the solve
        time (&end);
        dif = difftime (end,start);
 //       printf ("\nSolve took %.2lf seconds. \n", dif );
        
        double* voltage_array;
    
        // test whether voltages and gating variables are in correct ranges

        int ierr = VecGetArray(monodomainProblem.mCurrentVoltage, &voltage_array); 
        
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
            
            int num_nodes = monodomainProblem.mMesh.GetNumNodes();
            // test final voltages have returned to the resting potential
            for(int i=0; i<num_nodes; i++)
            {
                TS_ASSERT_DELTA(voltage_array[i], -84.5, 1);
            }
                        
        }
        
        VecRestoreArray(monodomainProblem.mCurrentVoltage, &voltage_array);      
        VecAssemblyBegin(monodomainProblem.mCurrentVoltage);
        VecAssemblyEnd(monodomainProblem.mCurrentVoltage);
        VecDestroy(monodomainProblem.mCurrentVoltage);
    }   
};
#endif //_TESTMONODOMAINDG0ASSEMBLERLONG_HPP_
