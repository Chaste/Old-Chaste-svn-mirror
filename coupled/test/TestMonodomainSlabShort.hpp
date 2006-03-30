#ifndef _TESTMONODOMAINSHORT_HPP_
#define _TESTMONODOMAINSHORT_HPP_

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
#include "SimpleDg0ParabolicAssembler.hpp"  
#include "MonodomainDg0Assembler.hpp"
#include "TrianglesMeshReader.hpp"
#include "ColumnDataWriter.hpp"
#include "ColumnDataReader.hpp"
#include "PropagationPropertiesCalculator.hpp"

#include "MonodomainPde.hpp"
#include "MockEulerIvpOdeSolver.hpp"
#include "FischerPde.hpp"

#include "PetscSetupAndFinalize.hpp"
#include "MonodomainProblem.hpp"
#include "AbstractLinearParabolicPde.hpp"
#include "AbstractMonodomainProblemStimulus.hpp"

class FaceStimulus3D: public AbstractMonodomainProblemStimulus<3>
{
    virtual void Apply(MonodomainPde<3> *pPde,
                       ConformingTetrahedralMesh<3,3> *pMesh)
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
};

class TestMonodomainSlab : public CxxTest::TestSuite 
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

    // Solve on a 3D 1mm by 1mm by 1mm mesh (space step = 0.1mm), stimulating 
    // the left face.
    // Should behave like the 1D case, extrapolated.
    // This version has a longer duration, and is disabled by default.  The
    // version below has the same duration as the tests in
    // TestMonodomainDg0Assembler.hpp

    void TestMonodomainDg03DWithFaceStimulus( void )
    {
        FaceStimulus3D face_stimulus_3D;
        
        MonodomainProblem<3> monodomainProblem;

        monodomainProblem.SetMeshFilename("mesh/test/data/3D_0_to_1mm_6000_elements");
        monodomainProblem.SetEndTime(5);   // 5 ms
        monodomainProblem.SetOutputDirectory("testoutput/MonoDg03dWithFaceStimulus");
        monodomainProblem.SetOutputFilenamePrefix("NewMonodomainLR91_3dWithFaceStimulus");
        monodomainProblem.SetStimulus(&face_stimulus_3D);

        clock_t start, stop;
        start=clock();
        monodomainProblem.Solve();
        stop=clock();
        double total_time=stop-start;
        
   //     TS_ASSERT_LESS_THAN(total_time/CLOCKS_PER_SEC, 18.0);
     // Change for 5 millisecond simulation time, ought to 2'32''=152 seconds
        TS_ASSERT_LESS_THAN(total_time/CLOCKS_PER_SEC, 160.0);
        
   //     std::cout<<"Total is "<<total_time/CLOCKS_PER_SEC<<"\n";
        
//        double* voltage_array;
//    
//        // test whether voltages and gating variables are in correct ranges
//
//        int ierr = VecGetArray(monodomainProblem.mCurrentVoltage, &voltage_array); 
//        
//        for(int global_index=monodomainProblem.mLo; global_index<monodomainProblem.mHi; global_index++)
//        {
//            // assuming LR model has Ena = 54.4 and Ek = -77
//            double Ena   =  54.4;
//            double Ek    = -77.0;
//            
//            TS_ASSERT_LESS_THAN_EQUALS(   voltage_array[global_index-monodomainProblem.mLo] , Ena +  30);
//            TS_ASSERT_LESS_THAN_EQUALS(  -voltage_array[global_index-monodomainProblem.mLo] + (Ek-30), 0);
//                
//            std::vector<double> odeVars = monodomainProblem.mMonodomainPde->GetOdeVarsAtNode(global_index);           
//            for(int j=0; j<8; j++)
//            {
//                // if not voltage or calcium ion conc, test whether between 0 and 1 
//                if((j!=4) && (j!=3))
//                {
//                    TS_ASSERT_LESS_THAN_EQUALS(  odeVars[j], 1.0);        
//                    TS_ASSERT_LESS_THAN_EQUALS( -odeVars[j], 0.0);        
//                }
//            }
//        }
//        
//        int num_procs;
//        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);
//
//        if (num_procs == 1)
//        {
//            /*
//             * Test the top right node against the right one in the 1D case, 
//             * comparing voltage, and then test all the nodes on the right hand 
//             * face of the cube against the top right one, comparing voltage.
//             */
//            bool need_initialisation = true;
//            double voltage;
//
//            need_initialisation = true;
//
//            // Test the RHF of the mesh
//            for (int i = 0; i < monodomainProblem.mMesh.GetNumNodes(); i++)
//            {
//                if (monodomainProblem.mMesh.GetNodeAt(i)->GetPoint()[0] == 0.1)
//                {
//                    // x = 0 is where the stimulus has been applied
//                    // x = 0.1cm is the other end of the mesh and where we want to 
//                    //       to test the value of the nodes
//                    
//                    if (need_initialisation)
//                    {
//                        voltage = voltage_array[i];
//                        need_initialisation = false;
//                    }
//                    else
//                    {
//                        TS_ASSERT_DELTA(voltage_array[i], voltage, 0.005);
//                       // std::cout << "y=" << monodomainProblem.mMesh.GetNodeAt(i)->GetPoint()[1] << std::endl;
//                    }
//                    
//                    // Check against 1d case
//                    TS_ASSERT_DELTA(voltage_array[i], 7.24231, 0.01);
//                }
//            }
//        }        
//        
//        VecRestoreArray(monodomainProblem.mCurrentVoltage, &voltage_array);      
//        VecAssemblyBegin(monodomainProblem.mCurrentVoltage);
//        VecAssemblyEnd(monodomainProblem.mCurrentVoltage);
//        VecDestroy(monodomainProblem.mCurrentVoltage);
    }   
};

#endif //_TESTMONODOMAINSHORT_HPP_
