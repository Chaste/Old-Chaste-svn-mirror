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
        
        // sets Luo Rudy system with initial conditions passed on
        double voltage = -9999; // This voltage will be ignored
        double m = 0.0017;
        double h = 0.9833;
        double j = 0.9895;
        double d = 0.003;
        double f = 1;
        double x = 0.0056;
        double caI = 0.0002; 
        double magnitudeOfStimulus = -80.0;  
        double durationOfStimulus  = 0.5 ;  // ms   
        
                  
        std::vector<double> initialConditions;
        initialConditions.push_back(h);
        initialConditions.push_back(j);
        initialConditions.push_back(m);
        initialConditions.push_back(caI);
        initialConditions.push_back(voltage);
        initialConditions.push_back(d);
        initialConditions.push_back(f);
        initialConditions.push_back(x);
        
        monodomain_pde.SetUniversalInitialConditions(initialConditions);
        
        InitialStimulus stimulus(magnitudeOfStimulus, durationOfStimulus);
           
        monodomain_pde.SetStimulusFunctionAtNode(0, &stimulus);
        
        // voltage that gets passed in solving ode
         voltage = -84.5;
 
   		// initial condition;   
		Vec currentVoltage;
		VecCreate(PETSC_COMM_WORLD, &currentVoltage);
		VecSetSizes(currentVoltage, PETSC_DECIDE, num_nodes);
		//VecSetType(initialCondition, VECSEQ);
		VecSetFromOptions(currentVoltage);
  
		double* currentVoltageArray;
		int ierr = VecGetArray(currentVoltage, &currentVoltageArray); 
        
		// initial voltage condition of a constant everywhere on the mesh
		
		currentVoltageArray[0] = -84.5;
		currentVoltageArray[1] = -84.5;
		
		VecRestoreArray(currentVoltage, &currentVoltageArray);      
		VecAssemblyBegin(currentVoltage);
		VecAssemblyEnd(currentVoltage);
		 
	    monodomain_pde.PrepareForAssembleSystem(currentVoltage);
        double value1 = monodomain_pde.ComputeNonlinearSourceTermAtNode(node0, voltage);
   
        initialConditions[4] = voltage;
        LuoRudyIModel1991OdeSystem Lr91OdeSystemStimulated(&stimulus);
                              
        OdeSolution SolutionNewStimulated = mySolver.Solve(&Lr91OdeSystemStimulated, start_time, start_time + big_time_step, small_time_step, initialConditions);  
        std::vector<double> solutionSetStimT_05 = SolutionNewStimulated.mSolutions[ SolutionNewStimulated.mSolutions.size()-1 ];
        
        double value2 = -(-80 + monodomain_pde.GetIIonic(solutionSetStimT_05));
        
        TS_ASSERT_DELTA(value1, value2, 0.000001);

        // shouldn't be different when called again as reset not yet been called
        value1 = monodomain_pde.ComputeNonlinearSourceTermAtNode(node0, voltage);
        TS_ASSERT_DELTA(value1, value2, 0.000001);
  
 
 
        InitialStimulus zeroStimulus(0, 0); 
        LuoRudyIModel1991OdeSystem Lr91OdeSystemNotStim(&zeroStimulus);

        OdeSolution SolutionNewNotStim = mySolver.Solve(&Lr91OdeSystemNotStim, start_time, start_time + big_time_step, small_time_step, initialConditions);  
        std::vector<double> solutionSetNoStimT_05 = SolutionNewNotStim.mSolutions[ SolutionNewNotStim.mSolutions.size()-1 ];
       
        value1 = monodomain_pde.ComputeNonlinearSourceTermAtNode(node1, voltage);
        value2 = -(0 + monodomain_pde.GetIIonic(solutionSetNoStimT_05));

        TS_ASSERT_DELTA(value1, value2, 0.000001);
 

 
 
        // Reset       
       	VecGetArray(currentVoltage, &currentVoltageArray); 
		currentVoltageArray[0] = solutionSetStimT_05[4];
		currentVoltageArray[1] = solutionSetNoStimT_05[4];
		
		VecRestoreArray(currentVoltage, &currentVoltageArray);      
		VecAssemblyBegin(currentVoltage);
		VecAssemblyEnd(currentVoltage);
		monodomain_pde.ResetAsUnsolvedOdeSystem();
        monodomain_pde.PrepareForAssembleSystem(currentVoltage);
              

        
        value1 = monodomain_pde.ComputeNonlinearSourceTermAtNode(node0, solutionSetStimT_05[4]);


        OdeSolution SolutionNewStimulatedT_1 = mySolver.Solve(&Lr91OdeSystemStimulated, start_time + big_time_step, 
                                                                start_time + 2*big_time_step, small_time_step,
                                                                solutionSetStimT_05);  
        std::vector<double> solutionSetStimT_1 = SolutionNewStimulatedT_1.mSolutions[ SolutionNewStimulatedT_1.mSolutions.size()-1 ];
        
        value2 = -(0 + monodomain_pde.GetIIonic(solutionSetStimT_1));
        
        TS_ASSERT_DELTA(value1, value2, 1e-10);
        
        




        OdeSolution SolutionNewNotStimT_1 = mySolver.Solve(&Lr91OdeSystemNotStim, start_time + big_time_step, 
                                                             start_time + 2*big_time_step, small_time_step, 
                                                             solutionSetNoStimT_05);
        
        std::vector<double> solutionSetNoStimT_1 = SolutionNewNotStimT_1.mSolutions[ SolutionNewNotStimT_1.mSolutions.size()-1 ];
       
        value1 = monodomain_pde.ComputeNonlinearSourceTermAtNode(node1, solutionSetNoStimT_05[4]);
        value2 = -(0 + monodomain_pde.GetIIonic(solutionSetNoStimT_1));

        TS_ASSERT_DELTA(value1, value2, 1e-10);

     

    }
    
     
};

#endif //_TESTMONODOMAINPDE_HPP_
