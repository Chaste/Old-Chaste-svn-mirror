#ifndef _TESTMONODOMAINPDE_HPP_
#define _TESTMONODOMAINPDE_HPP_

//#include <cmath>
#include <iostream>
//#include <fstream>
#include <vector>
#include "AbstractStimulusFunction.hpp"
//#include "InitialStimulus.hpp"
//#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "RungeKutta2IvpOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "AdamsBashforthIvpOdeSolver.hpp"
//#include "OdeSolution.hpp"
//#include "ColumnDataWriter.hpp"
//#include "LuoRudyIModel1991OdeSystem.hpp"

#include "MonodomainPde.hpp"


#include "EulerIvpOdeSolver.hpp"
#include "OdeSolution.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "petscvec.h"
 
// todo: test Fitzhugh Nagumo PDE 
 
#include <cxxtest/TestSuite.h>

class TestMonodomainPde : public CxxTest::TestSuite
{
    public:    
    void testMonodomainPde( void )
    {
        // Test for 2 nodes, check MonodomainPde correctly solves gating variable and ca_i concentration 
        // dynamics, comparing answers to LuoRudy data (chaste/data/Lr91Good.dat and Lr91NoStimGood.dat)
        // with no stimulus applied to node 1 and init stimulus applied to node 0. 
        // BigTimeStep = 1, SmallTimeStep = 0.01       
        // We really are solving extra ode (the one for voltage as its results are never used)
        int num_nodes=2; 
          
        Node<1> node0(0,true,0);
        Node<1> node1(1,true,0);
        
        double start_time = 0;  
        double big_time_step = 0.5;
        double small_time_step = 0.01; 
        EulerIvpOdeSolver mySolver;
        
        MonodomainPde<1> monodomain_pde(num_nodes, &mySolver, start_time, big_time_step, small_time_step);

        // Stimulus function to use at node 0. Node 1 is not stimulated.
        double magnitudeOfStimulus = -80.0;  
        double durationOfStimulus  = 0.5 ;  // ms   
        InitialStimulus stimulus(magnitudeOfStimulus, durationOfStimulus);
        monodomain_pde.SetStimulusFunctionAtNode(0, &stimulus);
        
        // voltage that gets passed in solving ode
        double voltage = -84.5;
 
   		// initial condition;   
		Vec currentVoltage;
		VecCreate(PETSC_COMM_WORLD, &currentVoltage);
		VecSetSizes(currentVoltage, PETSC_DECIDE, num_nodes);
		//VecSetType(initialCondition, VECSEQ);
		VecSetFromOptions(currentVoltage);
  
		double* currentVoltageArray;
		int ierr = VecGetArray(currentVoltage, &currentVoltageArray); 
        
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
   
        LuoRudyIModel1991OdeSystem ode_system_stimulated(&stimulus);
                              
        OdeSolution SolutionNewStimulated = mySolver.Solve(&ode_system_stimulated,
                                                           start_time,
                                                           start_time + big_time_step,
                                                           small_time_step,
                                                           ode_system_stimulated.mInitialConditions);  
        std::vector<double> solutionSetStimT_05 = SolutionNewStimulated.mSolutions[ SolutionNewStimulated.mSolutions.size()-1 ];
        
        double value2 = -(-80 + monodomain_pde.GetIIonic(solutionSetStimT_05));
        
        TS_ASSERT_DELTA(value1, value2, 0.000001);

        // shouldn't be different when called again as reset not yet been called
        value1 = monodomain_pde.ComputeNonlinearSourceTermAtNode(node0, voltage);
        TS_ASSERT_DELTA(value1, value2, 0.000001);
  
        InitialStimulus zeroStimulus(0, 0); 
        LuoRudyIModel1991OdeSystem ode_system_not_stim(&zeroStimulus);

        OdeSolution SolutionNewNotStim = mySolver.Solve(&ode_system_not_stim,
                                                        start_time,
                                                        start_time + big_time_step,
                                                        small_time_step,
                                                        ode_system_not_stim.mInitialConditions);  
        std::vector<double> solutionSetNoStimT_05 = SolutionNewNotStim.mSolutions[ SolutionNewNotStim.mSolutions.size()-1 ];
       
        value1 = monodomain_pde.ComputeNonlinearSourceTermAtNode(node1, voltage);
        value2 = -(0 + monodomain_pde.GetIIonic(solutionSetNoStimT_05));

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

        OdeSolution SolutionNewStimulatedT_1 = mySolver.Solve(&ode_system_stimulated, start_time + big_time_step, 
                                                                start_time + 2*big_time_step, small_time_step,
                                                                solutionSetStimT_05);  
        std::vector<double> solutionSetStimT_1 = SolutionNewStimulatedT_1.mSolutions[ SolutionNewStimulatedT_1.mSolutions.size()-1 ];
        
        value2 = -(0 + monodomain_pde.GetIIonic(solutionSetStimT_1));
        
        TS_ASSERT_DELTA(value1, value2, 1e-10);
        
        




        OdeSolution SolutionNewNotStimT_1 = mySolver.Solve(&ode_system_not_stim, start_time + big_time_step, 
                                                             start_time + 2*big_time_step, small_time_step, 
                                                             solutionSetNoStimT_05);
        
        std::vector<double> solutionSetNoStimT_1 = SolutionNewNotStimT_1.mSolutions[ SolutionNewNotStimT_1.mSolutions.size()-1 ];
       
        value1 = monodomain_pde.ComputeNonlinearSourceTermAtNode(node1, solutionSetNoStimT_05[4]);
        value2 = -(0 + monodomain_pde.GetIIonic(solutionSetNoStimT_1));

        TS_ASSERT_DELTA(value1, value2, 1e-10);

     
        VecDestroy(currentVoltage);
    }
    
     
};

#endif //_TESTMONODOMAINPDE_HPP_
