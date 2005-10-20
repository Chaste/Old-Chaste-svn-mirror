#ifndef _TESTMONODOMAINDG0ASSEMBLER_HPP_
#define _TESTMONODOMAINDG0ASSEMBLER_HPP_

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
};

class EdgeStimulus2D: public AbstractMonodomainProblemStimulus<2>
{
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
                       ConformingTetrahedralMesh<2,2> *pMesh)
    {
        static InitialStimulus stimulus(-600.0, 0.5);

        pPde->SetStimulusFunctionAtNode(mNode, &stimulus);
    }
};

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

class PointStimulus3D: public AbstractMonodomainProblemStimulus<3>
{
    virtual void Apply(MonodomainPde<3> *pPde,
                       ConformingTetrahedralMesh<3,3> *pMesh)
    {
        static InitialStimulus stimulus(-600.0, 0.5);

        pPde->SetStimulusFunctionAtNode(0, &stimulus);
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

    // Solve on a 1D string of cells, 1cm long.
	void TestMonodomainDg01D()
	{
        PointStimulus1D point_stimulus_1D;
        MonodomainProblem<1> monodomainProblem("mesh/test/data/1D_0_to_1_100_elements",
                                               30, // ms
                                               "testoutput/MonoDg01d",
                                               "NewMonodomainLR91_1d",
                                               &point_stimulus_1D);

        monodomainProblem.Solve();
        
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
               
            if (global_index==1)
            {
                TS_ASSERT_DELTA(voltage_array[global_index-monodomainProblem.mLo], 9.4435, 0.001);
            }
            if (global_index==3)
            {
                TS_ASSERT_DELTA(voltage_array[global_index-monodomainProblem.mLo], 9.4145, 0.001);
            }
            if (global_index==5)
            {
                TS_ASSERT_DELTA(voltage_array[global_index-monodomainProblem.mLo], 9.3630, 0.001);
            }
            if (global_index==7)
            {
                TS_ASSERT_DELTA(voltage_array[global_index-monodomainProblem.mLo], 9.2948, 0.001);
            }
            if (global_index==9)
            {
                TS_ASSERT_DELTA(voltage_array[global_index-monodomainProblem.mLo], 9.2135, 0.001);
            }
            if (global_index==75)
            {
                TS_ASSERT_DELTA(voltage_array[global_index-monodomainProblem.mLo], 5.4796, 0.001);
            }
        }





//        // Display the voltage at the right hand node of the mesh.
//        int num_procs;
//        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);
//
//        if (num_procs == 1)
//        {
//           for (int i = 0; i < monodomainProblem.mMesh.GetNumNodes(); i++)
//            {
//                if (monodomainProblem.mMesh.GetNodeAt(i)->GetPoint()[0] == 1.0)
//                {
//                    std::cout << "Voltage: " << voltage_array[i] << std::endl;
//                }
//            }
//        }





        VecRestoreArray(monodomainProblem.mCurrentVoltage, &voltage_array);      
        VecAssemblyBegin(monodomainProblem.mCurrentVoltage);
        VecAssemblyEnd(monodomainProblem.mCurrentVoltage);
        VecDestroy(monodomainProblem.mCurrentVoltage);

    }
    
    // Solve on a 2D 1cm by 1cm mesh, stimulating the left edge.
    // Should behave like the 1D case, extrapolated.
    // Not yet working.
    void xTestMonodomainDg02DWithEdgeStimulus( void )
    {   
        EdgeStimulus2D edge_stimulus_2D;
        
        MonodomainProblem<2> monodomainProblem("mesh/test/data/square_128_elements",
                                               0.1,
                                               "testoutput/MonoDg02dWithEdgeStimulus",
                                               "NewMonodomainLR91_2dWithEdgeStimulus",
                                               &edge_stimulus_2D);

        monodomainProblem.Solve();
        
        double* voltage_array;
        int ierr = VecGetArray(monodomainProblem.mCurrentVoltage, &voltage_array); 
    
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

            // Test the left hand side of the mesh, to check the voltage is
            // the 'same' all along it.
//            for (int i = 0; i < monodomainProblem.mMesh.GetNumNodes(); i++)
//            {
//                if (monodomainProblem.mMesh.GetNodeAt(i)->GetPoint()[0] == 0.0)
//                {
//                    if (need_initialisation)
//                    {
//                        voltage = voltage_array[i];
//                        need_initialisation = false;
//                    }
//                    else
//                    {
//                        TS_ASSERT_DELTA(voltage_array[i], voltage, 0.001);
//                        std::cout << "y=" << monodomainProblem.mMesh.GetNodeAt(i)->GetPoint()[1] << std::endl;
//                    }
//                }
//            }

            need_initialisation = true;

            // Test the RHS of the mesh
            for (int i = 0; i < monodomainProblem.mMesh.GetNumNodes(); i++)
            {
                if (monodomainProblem.mMesh.GetNodeAt(i)->GetPoint()[0] == 1.0)
                {
                    // x = 0 is where the stimulus has been applied
                    // x = 1 is the other end of the mesh and where we want to 
                    //       to test the value of the nodes
                    
                    if (need_initialisation)
                    {
                        voltage = voltage_array[i];
                        need_initialisation = false;
                    }
                    else
                    {
                        TS_ASSERT_DELTA(voltage_array[i], voltage, 0.001);
                       // std::cout << "y=" << monodomainProblem.mMesh.GetNodeAt(i)->GetPoint()[1] << std::endl;
                    }
                    
                    // Check against 1d case
                    TS_ASSERT_DELTA(voltage_array[i], 5.643, 0.001);
                }
            }
        }
        
        VecRestoreArray(monodomainProblem.mCurrentVoltage, &voltage_array);      
        VecAssemblyBegin(monodomainProblem.mCurrentVoltage);
        VecAssemblyEnd(monodomainProblem.mCurrentVoltage);
        VecDestroy(monodomainProblem.mCurrentVoltage);
    }   

    void TestMonodomainDg02DWithPointStimulus( void )
    {   
        PointStimulus2D point_stimulus_2D(4); // Central node
        
        MonodomainProblem<2> monodomainProblem("mesh/test/data/square_128_elements",
                                               0.1, 
                                               "testoutput/MonoDg02dWithPointStimulus",
                                               "NewMonodomainLR91_2dWithPointStimulus",
                                               &point_stimulus_2D);

        monodomainProblem.Solve();
        
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
             * Test that corners are equal, and centres of sides.
             */
            TS_ASSERT_DELTA(voltage_array[0], voltage_array[1], 0.0001);
            TS_ASSERT_DELTA(voltage_array[0], voltage_array[2], 0.0001);
            TS_ASSERT_DELTA(voltage_array[0], voltage_array[3], 0.0001);
            
            TS_ASSERT_DELTA(voltage_array[5], voltage_array[6], 0.0001);
            TS_ASSERT_DELTA(voltage_array[5], voltage_array[7], 0.0001);
            TS_ASSERT_DELTA(voltage_array[5], voltage_array[8], 0.0001);
        }
        
        VecRestoreArray(monodomainProblem.mCurrentVoltage, &voltage_array);      
        VecAssemblyBegin(monodomainProblem.mCurrentVoltage);
        VecAssemblyEnd(monodomainProblem.mCurrentVoltage);
        VecDestroy(monodomainProblem.mCurrentVoltage);
    }   

    void TestMonodomainDg03DWithFaceStimulus( void )
    {
        FaceStimulus3D face_stimulus_3D;
        
        MonodomainProblem<3> monodomainProblem("mesh/test/data/slab_395_elements",
                                               0.25, 
                                               "testoutput/MonoDg03dWithFaceStimulus",
                                               "NewMonodomainLR91_3dWithFaceStimulus",
                                               &face_stimulus_3D);
        ///todo This test now fails (not sensible ODE values) for mEndTime >= 0.26 

        monodomainProblem.Solve();
        
        double* voltage_array;
    
        // test whether voltages and gating variables are in correct ranges

        int ierr = VecGetArray(monodomainProblem.mCurrentVoltage, &voltage_array); 
        
        for(int global_index=monodomainProblem.mLo; global_index<monodomainProblem.mHi; global_index++)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77 and given magnitude of initial stim = -80
            double Ena   =  54.4;
            double Ek    = -77.0;
            double Istim = -80.0;
            
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
        VecRestoreArray(monodomainProblem.mCurrentVoltage, &voltage_array);      
        VecAssemblyBegin(monodomainProblem.mCurrentVoltage);
        VecAssemblyEnd(monodomainProblem.mCurrentVoltage);
        VecDestroy(monodomainProblem.mCurrentVoltage);
    }   

    void TestMonodomainDg03DWithPointStimulus( void )
    {
        PointStimulus3D point_stimulus_3D;
        
        MonodomainProblem<3> monodomainProblem("mesh/test/data/slab_395_elements",
                                               0.25, 
                                               "testoutput/MonoDg03dWithPointStimulus",
                                               "NewMonodomainLR91_3dWithPointStimulus",
                                               &point_stimulus_3D);
        ///todo This test now fails (not sensible ODE values) for mEndTime >= 0.26 

        monodomainProblem.Solve();
        
        double* voltage_array;
    
        // test whether voltages and gating variables are in correct ranges

        int ierr = VecGetArray(monodomainProblem.mCurrentVoltage, &voltage_array); 
        
        for(int global_index=monodomainProblem.mLo; global_index<monodomainProblem.mHi; global_index++)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77 and given magnitude of initial stim = -80
            double Ena   =  54.4;
            double Ek    = -77.0;
            double Istim = -80.0;
            
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
        VecRestoreArray(monodomainProblem.mCurrentVoltage, &voltage_array);      
        VecAssemblyBegin(monodomainProblem.mCurrentVoltage);
        VecAssemblyEnd(monodomainProblem.mCurrentVoltage);
        VecDestroy(monodomainProblem.mCurrentVoltage);
    }
    
    
    

};
#endif //_TESTMONODOMAINDG0ASSEMBLER_HPP_
