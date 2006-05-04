#ifndef _TESTMONODOMAINSLAB_HPP_
#define _TESTMONODOMAINSLAB_HPP_

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
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpStimulus, mTimeStep);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpZeroStimulus, mTimeStep);
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
    void longTestMonodomainDg03DWithFaceStimulus( void )
    {
        FaceStimulusCellFactory cell_factory;
        
        MonodomainProblem<3> monodomain_problem(&cell_factory);

        monodomain_problem.SetMeshFilename("mesh/test/data/3D_0_to_1mm_6000_elements");
        monodomain_problem.SetEndTime(60);   // 60 ms
        monodomain_problem.SetOutputDirectory("testoutput/MonoDg03dWithFaceStimulus");
        monodomain_problem.SetOutputFilenamePrefix("NewMonodomainLR91_3dWithFaceStimulus");
        monodomain_problem.Initialise();

        monodomain_problem.Solve();
        

        double* p_voltage_array;
        int lo, hi;
        monodomain_problem.GetVoltageArray(&p_voltage_array, lo, hi); 
    
        // test whether voltages and gating variables are in correct ranges
        for(int global_index=lo; global_index<hi; global_index++)
        {
            int local_index = global_index - lo;
            // assuming LR model has Ena = 54.4 and Ek = -77
            double Ena   =  54.4;
            double Ek    = -77.0;
            
            TS_ASSERT_LESS_THAN_EQUALS( p_voltage_array[local_index] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(-p_voltage_array[local_index] + (Ek-30), 0);
                
            std::vector<double> odeVars = monodomain_problem.GetMonodomainPde()->
                                           GetCardiacCell(global_index)->rGetStateVariables();
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
             * Test all the nodes on the right hand 
             * face of the cube against the top right one, comparing voltage.
             */
            bool need_initialisation = true;
            double voltage;

            need_initialisation = true;

            // Test the RHF of the mesh
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
                        TS_ASSERT_DELTA(p_voltage_array[i], voltage, 0.005);
                       // std::cout << "y=" << monodomain_problem.mMesh.GetNodeAt(i)->GetPoint()[1] << std::endl;
                    }
                }
            }
        }        
        monodomain_problem.RestoreVoltageArray( &p_voltage_array );
    }   


    // Solve on a 3D 1mm by 1mm by 1mm mesh (space step = 0.1mm), stimulating 
    // the left face.
    // Should behave like the 1D case, extrapolated.
    // See also TestMonodomainSlab.hpp
    void TestMonodomainDg03DWithFaceStimulus( void )
    {
        FaceStimulusCellFactory cell_factory;
        
        MonodomainProblem<3> monodomain_problem(&cell_factory);

        monodomain_problem.SetMeshFilename("mesh/test/data/3D_0_to_1mm_6000_elements");
        monodomain_problem.SetEndTime(4);   // 4 ms
        monodomain_problem.SetOutputDirectory("testoutput/MonoDg03dWithFaceStimulus");
        monodomain_problem.SetOutputFilenamePrefix("NewMonodomainLR91_3dWithFaceStimulus");
        monodomain_problem.Initialise();

        monodomain_problem.Solve();
        
        double* p_voltage_array;
        int lo, hi;
        monodomain_problem.GetVoltageArray(&p_voltage_array, lo, hi); 
    
        // test whether voltages and gating variables are in correct ranges
        for(int global_index=lo; global_index<hi; global_index++)
        {
            int local_index = global_index - lo;
            // assuming LR model has Ena = 54.4 and Ek = -77
            double Ena   =  54.4;
            double Ek    = -77.0;
            
            TS_ASSERT_LESS_THAN_EQUALS( p_voltage_array[local_index] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(-p_voltage_array[local_index] + (Ek-30), 0);
                
            std::vector<double> odeVars = monodomain_problem.GetMonodomainPde()->
                                           GetCardiacCell(global_index)->rGetStateVariables();
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
             * face of the cube against the top right one, comparing voltage.
             */
            bool need_initialisation = true;
            double voltage;

            need_initialisation = true;

            // Test the RHF of the mesh
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
                        TS_ASSERT_DELTA(p_voltage_array[i], voltage, 1);
                       // std::cout << "y=" << monodomain_problem.mMesh.GetNodeAt(i)->GetPoint()[1] << std::endl;
                    }
                    
                    // Check against 1d case - if the TestMonodomainDg01D test is run
                    // for 4ms the voltage at the end node is 21.8820
                    //TS_ASSERT_DELTA(p_voltage_array[i], 21.88, 0.5);
                    //
                    // This has now been changed as the initial conditions for Lr91 have been altered
                    TS_ASSERT_DELTA(p_voltage_array[i], 20.5, 0.5);
                    
                }
            }
        }        
        monodomain_problem.RestoreVoltageArray( &p_voltage_array );
    }   
};


#endif //_TESTMONODOMAINSLAB_HPP_
