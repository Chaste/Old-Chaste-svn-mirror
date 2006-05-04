#ifndef _TESTMONODOMAINPDE_HPP_
#define _TESTMONODOMAINPDE_HPP_

//#include <iostream>
#include <vector>

#include "InitialStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "MonodomainPde.hpp"
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
        

class TestMonodomainPde : public CxxTest::TestSuite
{
    public:    
    void testMonodomainPde( void )
    {

        int num_nodes=2; 
          
        Node<1> node0(0,true,0);
        Node<1> node1(1,true,0);
        
        double start_time = 0;  
        double big_time_step = 0.5;
        double small_time_step = 0.01; 
        
        AbstractIvpOdeSolver*   solver = new EulerIvpOdeSolver;
        InitialStimulus*     zero_stim = new InitialStimulus(0,0,0);

         // Stimulus function to use at node 0. Node 1 is not stimulated.
        double magnitudeOfStimulus = -80.0;  
        double durationOfStimulus  = 0.5 ;  // ms   

        AbstractStimulusFunction* stimulus = new InitialStimulus(magnitudeOfStimulus, durationOfStimulus);


        MyCardiacCellFactory cell_factory;

        MonodomainPde<1> monodomain_pde( &cell_factory, start_time, big_time_step );

        
        // voltage that gets passed in solving ode
        double initial_voltage = -83.853;
 
   		// initial condition;   
		Vec voltage;
		VecCreate(PETSC_COMM_WORLD, &voltage);
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


        double value1 = monodomain_pde.ComputeNonlinearSourceTermAtNode(node0, initial_voltage);
   
        LuoRudyIModel1991OdeSystem ode_system_stimulated(solver, stimulus, small_time_step);
                              
        OdeSolution SolutionNewStimulated = ode_system_stimulated.Compute(
                                                           start_time,
                                                           start_time + big_time_step);
        std::vector<double> solutionSetStimT_05 = SolutionNewStimulated.rGetSolutions()[ SolutionNewStimulated.rGetSolutions().size()-1 ];
        double value2 = -(-80 + ode_system_stimulated.GetIIonic());

        TS_ASSERT_DELTA(value1, value2, 0.000001);

        // shouldn't be different when called again as reset not yet been called
        value1 = monodomain_pde.ComputeNonlinearSourceTermAtNode(node0, initial_voltage);
        TS_ASSERT_DELTA(value1, value2, 0.000001);
  
        LuoRudyIModel1991OdeSystem ode_system_not_stim(solver, zero_stim, small_time_step);

        value1 = monodomain_pde.ComputeNonlinearSourceTermAtNode(node1, initial_voltage);

        OdeSolution SolutionNewNotStim = ode_system_not_stim.Compute(
                                                        start_time,
                                                        start_time + big_time_step);
        std::vector<double> solutionSetNoStimT_05 = SolutionNewNotStim.rGetSolutions()[ SolutionNewNotStim.rGetSolutions().size()-1 ];
        value2 = -(0 + ode_system_not_stim.GetIIonic());

        TS_ASSERT_DELTA(value1, value2, 0.000001);
 

        // Reset       
       	VecGetArray(voltage, &p_voltage_array); 
        

        if (lo<=0 && 0<hi)
        {
            int local_index = 0 - lo;
    		p_voltage_array[local_index] = solutionSetStimT_05[4];
        }
        if (lo<=1 && 1<hi)
        {
            int local_index = 1 - lo;
    		p_voltage_array[local_index] = solutionSetNoStimT_05[4];
        }
		
		VecRestoreArray(voltage, &p_voltage_array);      
		VecAssemblyBegin(voltage);
		VecAssemblyEnd(voltage);

		monodomain_pde.ResetAsUnsolvedOdeSystem();
        monodomain_pde.PrepareForAssembleSystem(voltage);
              
        value1 = monodomain_pde.ComputeNonlinearSourceTermAtNode(node0, solutionSetStimT_05[4]);

        std::vector<double> state_variables = solutionSetStimT_05;
        ode_system_stimulated.SetStateVariables(state_variables);
        OdeSolution SolutionNewStimulatedT_1 = ode_system_stimulated.Compute( start_time + big_time_step, start_time + 2*big_time_step );
        std::vector<double> solutionSetStimT_1 = SolutionNewStimulatedT_1.rGetSolutions()[ SolutionNewStimulatedT_1.rGetSolutions().size()-1 ];
        value2 = -(0 + ode_system_stimulated.GetIIonic());
                
        TS_ASSERT_DELTA(value1, value2, 1e-10);
        
        state_variables = solutionSetNoStimT_05;
        ode_system_not_stim.SetStateVariables(state_variables);
        OdeSolution SolutionNewNotStimT_1 = ode_system_not_stim.Compute( start_time + big_time_step, start_time + 2*big_time_step );
        std::vector<double> solutionSetNoStimT_1 = SolutionNewNotStimT_1.rGetSolutions()[ SolutionNewNotStimT_1.rGetSolutions().size()-1 ];
       
        value1 = monodomain_pde.ComputeNonlinearSourceTermAtNode(node1, solutionSetNoStimT_05[4]);
        value2 = -(0 + ode_system_not_stim.GetIIonic());
        
        
        TS_ASSERT_DELTA(value1, value2, 1e-10);

     
        VecDestroy(voltage);
        delete zero_stim;
        delete stimulus;
        delete solver;        
    }
    
    
    void testMonodomainPdeGetCardiacCell( void )
    {
        int num_nodes = 2;
        MyCardiacCellFactory cell_factory;
        MonodomainPde<1> monodomain_pde( &cell_factory, 0, 0.1 );
        
        // initial condition;   
        Vec voltage;
        VecCreate(PETSC_COMM_WORLD, &voltage);
        VecSetSizes(voltage, PETSC_DECIDE, num_nodes);
        //VecSetType(initialCondition, VECSEQ);
        VecSetFromOptions(voltage);
  
        double* p_voltage_array;
        VecGetArray(voltage, &p_voltage_array); 
        
        int lo, hi;
        VecGetOwnershipRange(voltage,&lo,&hi);

        if(lo<=0 && 0<hi)
        {
            AbstractCardiacCell* cell = monodomain_pde.GetCardiacCell(0);
            TS_ASSERT_DELTA(cell->GetStimulus(0.001),-80,1e-10);
        }
    
        if(lo<=1 && 1<hi)
        {
            AbstractCardiacCell* cell = monodomain_pde.GetCardiacCell(1);
            TS_ASSERT_DELTA(cell->GetStimulus(0.001),0,1e-10);
        }
        
        VecDestroy(voltage);
    }
};



#endif //_TESTMONODOMAINPDE_HPP_
