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

class PointStimulus1D: public AbstractMonodomainProblemStimulus<1>
{
public:
    virtual void Apply(MonodomainPde<1> *pPde)
    {
        static InitialStimulus stimulus(-600.0, 0.5);
        pPde->SetStimulusFunctionAtNode(0, &stimulus);
    }
};

class EdgeStimulus2D: public AbstractMonodomainProblemStimulus<2>
{
    virtual void Apply(MonodomainPde<2> *pPde)
    {
        static InitialStimulus stimulus(-600.0, 0.5);
        pPde->SetStimulusFunctionAtNode(0, &stimulus);
        pPde->SetStimulusFunctionAtNode(3, &stimulus);
        pPde->SetStimulusFunctionAtNode(5, &stimulus);
        pPde->SetStimulusFunctionAtNode(14, &stimulus);
        pPde->SetStimulusFunctionAtNode(17, &stimulus);
        pPde->SetStimulusFunctionAtNode(43, &stimulus);
        pPde->SetStimulusFunctionAtNode(52, &stimulus);
        pPde->SetStimulusFunctionAtNode(53, &stimulus);
        pPde->SetStimulusFunctionAtNode(63, &stimulus);
    }
};

class PointStimulus2D: public AbstractMonodomainProblemStimulus<2>
{
    virtual void Apply(MonodomainPde<2> *pPde)
    {
        static InitialStimulus stimulus(-600.0, 0.5);
        pPde->SetStimulusFunctionAtNode(0, &stimulus);
    }
};

class FaceStimulus3D: public AbstractMonodomainProblemStimulus<3>
{
    virtual void Apply(MonodomainPde<3> *pPde)
    {
        static InitialStimulus stimulus(-600.0, 0.5);
        pPde->SetStimulusFunctionAtNode(0, &stimulus);
        pPde->SetStimulusFunctionAtNode(3, &stimulus);
        pPde->SetStimulusFunctionAtNode(4, &stimulus);
        pPde->SetStimulusFunctionAtNode(7, &stimulus);
        pPde->SetStimulusFunctionAtNode(8, &stimulus);
        pPde->SetStimulusFunctionAtNode(11, &stimulus);
        pPde->SetStimulusFunctionAtNode(17, &stimulus);
        pPde->SetStimulusFunctionAtNode(18, &stimulus);
        pPde->SetStimulusFunctionAtNode(23, &stimulus);
        pPde->SetStimulusFunctionAtNode(32, &stimulus);
        pPde->SetStimulusFunctionAtNode(33, &stimulus);
        pPde->SetStimulusFunctionAtNode(34, &stimulus);
        pPde->SetStimulusFunctionAtNode(41, &stimulus);
        pPde->SetStimulusFunctionAtNode(45, &stimulus);
        pPde->SetStimulusFunctionAtNode(46, &stimulus);
        pPde->SetStimulusFunctionAtNode(52, &stimulus);
        pPde->SetStimulusFunctionAtNode(56, &stimulus);
        pPde->SetStimulusFunctionAtNode(57, &stimulus);
        pPde->SetStimulusFunctionAtNode(58, &stimulus);
        pPde->SetStimulusFunctionAtNode(78, &stimulus);
        pPde->SetStimulusFunctionAtNode(80, &stimulus);
        pPde->SetStimulusFunctionAtNode(83, &stimulus);
        pPde->SetStimulusFunctionAtNode(86, &stimulus);
        pPde->SetStimulusFunctionAtNode(87, &stimulus);
        pPde->SetStimulusFunctionAtNode(88, &stimulus);
        pPde->SetStimulusFunctionAtNode(106, &stimulus);
        pPde->SetStimulusFunctionAtNode(113, &stimulus);
    }
};

class PointStimulus3D: public AbstractMonodomainProblemStimulus<3>
{
    virtual void Apply(MonodomainPde<3> *pPde)
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
	void TestMonodomainDg01D()
	{
        PointStimulus1D point_stimulus_1D;
        MonodomainProblem<1> monodomainProblem("mesh/test/data/1D_0_to_1_10_elements",
                                               30, 
                                               "testoutput/MonoDg01d",
                                               "NewMonodomainLR91_1d",
                                               &point_stimulus_1D);

        monodomainProblem.Solve();
        
        double* currentVoltageArray;
    
        // test whether voltages and gating variables are in correct ranges

        int ierr = VecGetArray(monodomainProblem.mCurrentVoltage, &currentVoltageArray); 
        
        for(int global_index=monodomainProblem.mLo; global_index<monodomainProblem.mHi; global_index++)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77 and given magnitude of initial stim = -80
            double Ena   =  54.4;
            double Ek    = -77.0;
            double Istim = -80.0;
            
            TS_ASSERT_LESS_THAN_EQUALS(   currentVoltageArray[global_index-monodomainProblem.mLo] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(  -currentVoltageArray[global_index-monodomainProblem.mLo] + (Ek-30), 0);
                
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
                TS_ASSERT_DELTA(currentVoltageArray[global_index-monodomainProblem.mLo], 7.4940, 0.001);
            }
            if (global_index==3)
            {
                TS_ASSERT_DELTA(currentVoltageArray[global_index-monodomainProblem.mLo], 0.3840, 0.001);
            }
            if (global_index==5)
            {
                TS_ASSERT_DELTA(currentVoltageArray[global_index-monodomainProblem.mLo], -52.8798, 0.001);
            }
            if (global_index==7)
            {
                TS_ASSERT_DELTA(currentVoltageArray[global_index-monodomainProblem.mLo], -84.9598, 0.001);
            }
            if (global_index==9)
            {
                TS_ASSERT_DELTA(currentVoltageArray[global_index-monodomainProblem.mLo], -84.5556, 0.001);
            }
            
        }
        VecRestoreArray(monodomainProblem.mCurrentVoltage, &currentVoltageArray);      
        VecAssemblyBegin(monodomainProblem.mCurrentVoltage);
        VecAssemblyEnd(monodomainProblem.mCurrentVoltage);
        VecDestroy(monodomainProblem.mCurrentVoltage);

    }
    
    void TestMonodomainDg02DWithEdgeStimulus( void )
    {   
        EdgeStimulus2D edge_stimulus_2D;
        
        MonodomainProblem<2> monodomainProblem("mesh/test/data/square_128_elements",
                                               0.1, 
                                               "testoutput/MonoDg02dWithEdgeStimulus",
                                               "NewMonodomainLR91_2dWithEdgeStimulus",
                                               &edge_stimulus_2D);

        monodomainProblem.Solve();
        
        double* currentVoltageArray;
    
        // test whether voltages and gating variables are in correct ranges

        int ierr = VecGetArray(monodomainProblem.mCurrentVoltage, &currentVoltageArray); 
        
        for(int global_index=monodomainProblem.mLo; global_index<monodomainProblem.mHi; global_index++)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77 and given magnitude of initial stim = -80
            double Ena   =  54.4;
            double Ek    = -77.0;
            double Istim = -80.0;
            
            TS_ASSERT_LESS_THAN_EQUALS(   currentVoltageArray[global_index-monodomainProblem.mLo] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(  -currentVoltageArray[global_index-monodomainProblem.mLo] + (Ek-30), 0);
                
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
             * comparing voltage
             */

// TO BE DONE!
//            TS_ASSERT_DELTA(currentVoltageArray[1], <some value>, 0.001);
             
            /*
             * Test all the nodes on the right hand side of the square against 
             * the top right one, comparing voltage
             */
            TS_ASSERT_DELTA(currentVoltageArray[2], currentVoltageArray[1], 0.001);
//            TS_ASSERT_DELTA(currentVoltageArray[6], currentVoltageArray[1], 0.001);
//            TS_ASSERT_DELTA(currentVoltageArray[21], currentVoltageArray[1], 0.001);
//            TS_ASSERT_DELTA(currentVoltageArray[23], currentVoltageArray[1], 0.001);
//            TS_ASSERT_DELTA(currentVoltageArray[73], currentVoltageArray[1], 0.001);
//            TS_ASSERT_DELTA(currentVoltageArray[76], currentVoltageArray[1], 0.001);
//            TS_ASSERT_DELTA(currentVoltageArray[77], currentVoltageArray[1], 0.001);
//            TS_ASSERT_DELTA(currentVoltageArray[79], currentVoltageArray[1], 0.001);
        }
        
        VecRestoreArray(monodomainProblem.mCurrentVoltage, &currentVoltageArray);      
        VecAssemblyBegin(monodomainProblem.mCurrentVoltage);
        VecAssemblyEnd(monodomainProblem.mCurrentVoltage);
        VecDestroy(monodomainProblem.mCurrentVoltage);
    }   

    void TestMonodomainDg02DWithPointStimulus( void )
    {   
        PointStimulus2D point_stimulus_2D;
        
        MonodomainProblem<2> monodomainProblem("mesh/test/data/square_128_elements",
                                               0.1, 
                                               "testoutput/MonoDg02dWithPointStimulus",
                                               "NewMonodomainLR91_2dWithPointStimulus",
                                               &point_stimulus_2D);

        monodomainProblem.Solve();
        
        double* currentVoltageArray;
    
        // test whether voltages and gating variables are in correct ranges

        int ierr = VecGetArray(monodomainProblem.mCurrentVoltage, &currentVoltageArray); 
        
        for(int global_index=monodomainProblem.mLo; global_index<monodomainProblem.mHi; global_index++)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77 and given magnitude of initial stim = -80
            double Ena   =  54.4;
            double Ek    = -77.0;
            double Istim = -80.0;
            
            TS_ASSERT_LESS_THAN_EQUALS(   currentVoltageArray[global_index-monodomainProblem.mLo] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(  -currentVoltageArray[global_index-monodomainProblem.mLo] + (Ek-30), 0);
                
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
        VecRestoreArray(monodomainProblem.mCurrentVoltage, &currentVoltageArray);      
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
        
        double* currentVoltageArray;
    
        // test whether voltages and gating variables are in correct ranges

        int ierr = VecGetArray(monodomainProblem.mCurrentVoltage, &currentVoltageArray); 
        
        for(int global_index=monodomainProblem.mLo; global_index<monodomainProblem.mHi; global_index++)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77 and given magnitude of initial stim = -80
            double Ena   =  54.4;
            double Ek    = -77.0;
            double Istim = -80.0;
            
            TS_ASSERT_LESS_THAN_EQUALS(   currentVoltageArray[global_index-monodomainProblem.mLo] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(  -currentVoltageArray[global_index-monodomainProblem.mLo] + (Ek-30), 0);
                
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
        VecRestoreArray(monodomainProblem.mCurrentVoltage, &currentVoltageArray);      
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
        
        double* currentVoltageArray;
    
        // test whether voltages and gating variables are in correct ranges

        int ierr = VecGetArray(monodomainProblem.mCurrentVoltage, &currentVoltageArray); 
        
        for(int global_index=monodomainProblem.mLo; global_index<monodomainProblem.mHi; global_index++)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77 and given magnitude of initial stim = -80
            double Ena   =  54.4;
            double Ek    = -77.0;
            double Istim = -80.0;
            
            TS_ASSERT_LESS_THAN_EQUALS(   currentVoltageArray[global_index-monodomainProblem.mLo] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(  -currentVoltageArray[global_index-monodomainProblem.mLo] + (Ek-30), 0);
                
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
        VecRestoreArray(monodomainProblem.mCurrentVoltage, &currentVoltageArray);      
        VecAssemblyBegin(monodomainProblem.mCurrentVoltage);
        VecAssemblyEnd(monodomainProblem.mCurrentVoltage);
        VecDestroy(monodomainProblem.mCurrentVoltage);
    }   
};
#endif //_TESTMONODOMAINDG0ASSEMBLER_HPP_
