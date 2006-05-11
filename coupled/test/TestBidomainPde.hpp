#ifndef TESTBIDOMAINPDE_HPP_
#define TESTBIDOMAINPDE_HPP_


//#include <iostream>
#include <vector>

#include "InitialStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "MonodomainPde.hpp"
#include "BidomainPde.hpp"
#include "OdeSolution.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include <petsc.h>
#include <cxxtest/TestSuite.h>



class MyCardiacCellFactory : public AbstractCardiacCellFactory<1>
{
private:
   InitialStimulus* mpStimulus;

public:
    
    MyCardiacCellFactory() : AbstractCardiacCellFactory<1>(0.01)
    {
        mpStimulus = new InitialStimulus(-80.0, 0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(int node)
    {                    
        if(node==0)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpStimulus, mTimeStep);
        }
        else if(node==1)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mpZeroStimulus, mTimeStep);
        }
        else
        {
            assert(0);
        }

    }
    
    ~MyCardiacCellFactory(void)
    {
        delete mpStimulus;
    }
    
    int GetNumberOfNodes()
    {
        return 2;
    }
};



class TestBidomainPde : public CxxTest::TestSuite
{
    public:    
    void testBidomainPde( void )
    {   
        double start_time = 0;  
        double big_time_step = 0.5;        
        MyCardiacCellFactory cell_factory;

        MonodomainPde<1> monodomain_pde( &cell_factory, start_time, big_time_step );        
        BidomainPde<1>   bidomain_pde( &cell_factory, start_time, big_time_step );   
        
        // voltage that gets passed in solving ode
        double initial_voltage = -83.853;
 
        // initial condition;   
        Vec voltage;
        VecCreate(PETSC_COMM_WORLD, &voltage);
        
        int num_nodes = 2;
        
        VecSetSizes(voltage, PETSC_DECIDE, num_nodes);
        //VecSetType(initialCondition, VECSEQ);
        VecSetFromOptions(voltage);
  
        double* p_voltage_array;
        VecGetArray(voltage, &p_voltage_array); 
        
        int lo, hi;
        VecGetOwnershipRange(voltage,&lo,&hi);
        
        for (int global_index = 0; global_index < num_nodes; global_index++ )
        {
            int local_index = global_index - lo;
            // initial voltage condition of a constant everywhere on the mesh
            if (lo<=global_index && global_index<hi)
            {
                p_voltage_array[local_index] = initial_voltage;
            }
        }
        VecRestoreArray(voltage, &p_voltage_array);      
        VecAssemblyBegin(voltage);
        VecAssemblyEnd(voltage);
         
       
        monodomain_pde.PrepareForAssembleSystem(voltage);
        bidomain_pde.PrepareForAssembleSystem(voltage);         

        // Check that both the monodomain and bidomain PDE have the same ionic cache

        for (int global_index=lo; global_index < hi; global_index++)
        {
            TS_ASSERT_EQUALS(monodomain_pde.GetIionicCacheReplicated()[global_index], bidomain_pde.GetIionicCacheReplicated()[global_index]);
        }

        // Check that the bidomain PDE has the right intracellular stimulus at node 0 and 1

        TS_ASSERT_EQUALS(bidomain_pde.GetIntracellularStimulusCacheReplicated()[0], -80);
        TS_ASSERT_EQUALS(bidomain_pde.GetIntracellularStimulusCacheReplicated()[1], 0);

        // TODO: check that the bidomain PDE has the right extracellular stimulus at node 0 and 1
        // NOTE: this will involve creating a new "MyCardiacCellFactory" that has both 
        //       an intracellular and an extracellular stimulus. The monodomain PDE
        //       would use the "MyCardiacCellFactory", while the bidomain PDE would use
        //       the new one
    }
};        

#endif /*TESTBIDOMAINPDE_HPP_*/
