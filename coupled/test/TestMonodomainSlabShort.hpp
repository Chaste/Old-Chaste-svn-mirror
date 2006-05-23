#ifndef _TESTMONODOMAINSHORT_HPP_
#define _TESTMONODOMAINSHORT_HPP_

// Element.hpp includes the Boost ublas objects - these need to
// be included early...  We think.  We're not that sure.
#include "Element.hpp"

#include <cxxtest/TestSuite.h>
#include "petscvec.h"
#include <vector>
//#include <iostream>

#include "ConformingTetrahedralMesh.cpp"
#include "PetscSetupAndFinalize.hpp"
#include "MonodomainProblem.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"


class FaceStimulusCellFactory : public AbstractCardiacCellFactory<3>
{
private:
    InitialStimulus *mpStimulus;
public:
    FaceStimulusCellFactory() : AbstractCardiacCellFactory<3>(0.01)
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
    
    ~FaceStimulusCellFactory(void)
    {
        delete mpStimulus;
    }
};



class TestMonodomainSlab : public CxxTest::TestSuite 
{   
public:

    // Solve on a 3D 1mm by 1mm by 1mm mesh (space step = 0.1mm), stimulating 
    // the left face.
    // Should behave like the 1D case, extrapolated.
    // This version has a longer duration, and is disabled by default.  The
    // version below has the same duration as the tests in
    // TestMonodomainDg0Assembler.hpp

    void TestMonodomainDg03DWithFaceStimulus( void )
    {
        FaceStimulusCellFactory cell_factory;
        
        MonodomainProblem<3> monodomain_problem(&cell_factory);

        monodomain_problem.SetMeshFilename("mesh/test/data/3D_0_to_1mm_6000_elements");
        monodomain_problem.SetEndTime(5);   // 5 ms
        monodomain_problem.SetOutputDirectory("testoutput/MonoDg03dWithFaceStimulus");
        monodomain_problem.SetOutputFilenamePrefix("NewMonodomainLR91_3dWithFaceStimulus");
        monodomain_problem.Initialise();

        clock_t start, stop;
        start=clock();
        monodomain_problem.Solve();
        stop=clock();
        double total_time=stop-start;
        
        //TS_ASSERT_LESS_THAN(total_time/CLOCKS_PER_SEC, 18.0);
        // Change for 5 millisecond simulation time, ought to 2'32''=152 seconds
        TS_ASSERT_LESS_THAN(total_time/CLOCKS_PER_SEC, 160.0);
        
//        std::cout<<"Total is "<<total_time/CLOCKS_PER_SEC<<"\n";
        
//        double* voltage_array;
//    
//        // test whether voltages and gating variables are in correct ranges
//
//        int ierr = VecGetArray(monodomain_problem.mCurrentVoltage, &voltage_array); 
//        
//        for(int global_index=monodomain_problem.mLo; global_index<monodomain_problem.mHi; global_index++)
//        {
//            // assuming LR model has Ena = 54.4 and Ek = -77
//            double Ena   =  54.4;
//            double Ek    = -77.0;
//            
//            TS_ASSERT_LESS_THAN_EQUALS(   voltage_array[global_index-monodomain_problem.mLo] , Ena +  30);
//            TS_ASSERT_LESS_THAN_EQUALS(  -voltage_array[global_index-monodomain_problem.mLo] + (Ek-30), 0);
//                
//            std::vector<double> odeVars = monodomain_problem.mMonodomainPde->GetCardiacCell(global_index)->GetStateVariables();
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
//            for (int i = 0; i < monodomain_problem.mMesh.GetNumNodes(); i++)
//            {
//                if (monodomain_problem.mMesh.GetNodeAt(i)->GetPoint()[0] == 0.1)
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
//                       // std::cout << "y=" << monodomain_problem.mMesh.GetNodeAt(i)->GetPoint()[1] << std::endl;
//                    }
//                    
//                    // Check against 1d case
//                    TS_ASSERT_DELTA(voltage_array[i], 7.24231, 0.01);
//                }
//            }
//        }        
//        
//        VecRestoreArray(monodomain_problem.mCurrentVoltage, &voltage_array);      
//        VecAssemblyBegin(monodomain_problem.mCurrentVoltage);
//        VecAssemblyEnd(monodomain_problem.mCurrentVoltage);
//        VecDestroy(monodomain_problem.mCurrentVoltage);
    }   
};

#endif //_TESTMONODOMAINSHORT_HPP_
