#ifndef _TESTMONODOMAINPDE_HPP_
#define _TESTMONODOMAINPDE_HPP_

//#include <iostream>
#include <vector>
#include "InitialStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"

#include "LuoRudyIModel1991OdeSystem.hpp"
#include "MonodomainPdeIteration7.hpp"

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
    
    MyCardiacCellFactory()
    {
        mTimeStep = 0.01;
        mpSolver = new EulerIvpOdeSolver;
        mpZeroStimulus = new InitialStimulus(0,0,0);
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
        delete mpSolver;
        delete mpZeroStimulus;
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

        MonodomainPdeIteration7<1> monodomain_pde( &cell_factory, start_time, big_time_step );

        
        // voltage that gets passed in solving ode
        double voltage = -84.5;
 
   		// initial condition;   
		Vec currentVoltage;
		VecCreate(PETSC_COMM_WORLD, &currentVoltage);
		VecSetSizes(currentVoltage, PETSC_DECIDE, num_nodes);
		//VecSetType(initialCondition, VECSEQ);
		VecSetFromOptions(currentVoltage);
  
		double* currentVoltageArray;
		VecGetArray(currentVoltage, &currentVoltageArray); 
        
        int lo, hi;
        VecGetOwnershipRange(currentVoltage,&lo,&hi);
        
		// initial voltage condition of a constant everywhere on the mesh
		if (lo<=0 && 0<hi)
        {
            currentVoltageArray[0-lo] = voltage;
        }
        if (lo<=1 && 1<hi)
        {
		  currentVoltageArray[1-lo] = voltage;
        }
		
		VecRestoreArray(currentVoltage, &currentVoltageArray);      
		VecAssemblyBegin(currentVoltage);
		VecAssemblyEnd(currentVoltage);
		 
	    monodomain_pde.PrepareForAssembleSystem(currentVoltage);


        double value1 = monodomain_pde.ComputeNonlinearSourceTermAtNode(node0, voltage);
   
        LuoRudyIModel1991OdeSystem ode_system_stimulated(solver, stimulus, small_time_step);
                              
        OdeSolution SolutionNewStimulated = ode_system_stimulated.Compute(
                                                           start_time,
                                                           start_time + big_time_step);
        std::vector<double> solutionSetStimT_05 = SolutionNewStimulated.mSolutions[ SolutionNewStimulated.mSolutions.size()-1 ];
        double value2 = -(-80 + ode_system_stimulated.GetIIonic());

        TS_ASSERT_DELTA(value1, value2, 0.000001);

        // shouldn't be different when called again as reset not yet been called
        value1 = monodomain_pde.ComputeNonlinearSourceTermAtNode(node0, voltage);
        TS_ASSERT_DELTA(value1, value2, 0.000001);
  
        LuoRudyIModel1991OdeSystem ode_system_not_stim(solver, zero_stim, small_time_step);

        value1 = monodomain_pde.ComputeNonlinearSourceTermAtNode(node1, voltage);

        OdeSolution SolutionNewNotStim = ode_system_not_stim.Compute(
                                                        start_time,
                                                        start_time + big_time_step);
        std::vector<double> solutionSetNoStimT_05 = SolutionNewNotStim.mSolutions[ SolutionNewNotStim.mSolutions.size()-1 ];
        value2 = -(0 + ode_system_not_stim.GetIIonic());

        TS_ASSERT_DELTA(value1, value2, 0.000001);
 

        // Reset       
       	VecGetArray(currentVoltage, &currentVoltageArray); 
        
        if (lo<=0 && 0<hi)
        {
    		currentVoltageArray[0-lo] = solutionSetStimT_05[4];
        }
        if (lo<=1 && 1<hi)
        {
    		currentVoltageArray[1-lo] = solutionSetNoStimT_05[4];
        }
		
		VecRestoreArray(currentVoltage, &currentVoltageArray);      
		VecAssemblyBegin(currentVoltage);
		VecAssemblyEnd(currentVoltage);

		monodomain_pde.ResetAsUnsolvedOdeSystem();
        monodomain_pde.PrepareForAssembleSystem(currentVoltage);
              
        value1 = monodomain_pde.ComputeNonlinearSourceTermAtNode(node0, solutionSetStimT_05[4]);

        std::vector<double> state_variables = solutionSetStimT_05;
        ode_system_stimulated.SetStateVariables(state_variables);
        OdeSolution SolutionNewStimulatedT_1 = ode_system_stimulated.Compute( start_time + big_time_step, start_time + 2*big_time_step );
        std::vector<double> solutionSetStimT_1 = SolutionNewStimulatedT_1.mSolutions[ SolutionNewStimulatedT_1.mSolutions.size()-1 ];
        value2 = -(0 + ode_system_stimulated.GetIIonic());
                
        TS_ASSERT_DELTA(value1, value2, 1e-10);
        
        state_variables = solutionSetNoStimT_05;
        ode_system_not_stim.SetStateVariables(state_variables);
        OdeSolution SolutionNewNotStimT_1 = ode_system_not_stim.Compute( start_time + big_time_step, start_time + 2*big_time_step );
        std::vector<double> solutionSetNoStimT_1 = SolutionNewNotStimT_1.mSolutions[ SolutionNewNotStimT_1.mSolutions.size()-1 ];
       
        value1 = monodomain_pde.ComputeNonlinearSourceTermAtNode(node1, solutionSetNoStimT_05[4]);
        value2 = -(0 + ode_system_not_stim.GetIIonic());
        
        
        TS_ASSERT_DELTA(value1, value2, 1e-10);

     
        VecDestroy(currentVoltage);
        delete zero_stim;
        delete stimulus;
        delete solver;        
    }
};

#endif //_TESTMONODOMAINPDE_HPP_
