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
        AbstractIvpOdeSolver *pMySolver = new EulerIvpOdeSolver();
        
        MonodomainPde<1> monodomain_pde(num_nodes, pMySolver, start_time, big_time_step, small_time_step);
        
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
        
        AbstractStimulusFunction *pStimulus = new InitialStimulus(magnitudeOfStimulus, durationOfStimulus);
           
        monodomain_pde.SetStimulusFunctionAtNode(0, pStimulus);
        
        // voltage that gets passed in solving ode
        voltage = -84.5;
        
        
/*      double i_stim = 0.0; // Set i_stim to zero by default
        
        double i_ionic;
        double i_total;        
                  
        //  Comparing data against data for known solution: t = 0.5
        double v = -4.4120e+01 ;
        m = 2.4052e-01;
        h = 9.6746e-01 ; 
        j = 9.8200e-01; 
        d = 3.6459e-03;
        f = 9.9997e-01;
        x = 5.6535e-03;
        caI = 1.9913e-04;
                  
        std::vector<double> solutionSetStimT_05;
        solutionSetStimT_05.push_back(h);
        solutionSetStimT_05.push_back(j);
        solutionSetStimT_05.push_back(m);
        solutionSetStimT_05.push_back(caI);
        solutionSetStimT_05.push_back(v);
        solutionSetStimT_05.push_back(d);
        solutionSetStimT_05.push_back(f);
        solutionSetStimT_05.push_back(x);
        
        magnitudeOfStimulus = -80; //at time t=1, stimulus is zero
        i_stim = magnitudeOfStimulus;
        
         // voltage that gets passed in solving ode
        voltage = -84.5;
        
        i_ionic = monodomain_pde.GetIIonic(solutionSetStimT_05, voltage);
        i_total = i_stim + i_ionic;
  */ 
  
        double value1 = monodomain_pde.ComputeNonlinearSourceTermAtNode(node0, voltage);
   
        initialConditions[4] = voltage;
        AbstractOdeSystem *pLr91OdeSystemStimulated = new LuoRudyIModel1991OdeSystem(pStimulus);
                              
        OdeSolution SolutionNewStimulated = pMySolver->Solve(pLr91OdeSystemStimulated, start_time, start_time + big_time_step, small_time_step, initialConditions);  
        std::vector<double> solutionSetStimT_05 = SolutionNewStimulated.mSolutions[ SolutionNewStimulated.mSolutions.size()-1 ];
        
        double value2 = -(-80 + monodomain_pde.GetIIonic(solutionSetStimT_05, voltage));
        
        TS_ASSERT_DELTA(value1, value2, 0.000001);

        // shouldn't be different when called again as reset not yet been called
        value1 = monodomain_pde.ComputeNonlinearSourceTermAtNode(node0, voltage);
        TS_ASSERT_DELTA(value1, value2, 0.000001);
  
 
 
        AbstractStimulusFunction* pZeroStimulus = new InitialStimulus(0, 0); 
        AbstractOdeSystem *pLr91OdeSystemNotStim = new LuoRudyIModel1991OdeSystem(pZeroStimulus);

        OdeSolution SolutionNewNotStim = pMySolver->Solve(pLr91OdeSystemNotStim, start_time, start_time + big_time_step, small_time_step, initialConditions);  
        std::vector<double> solutionSetNoStimT_05 = SolutionNewNotStim.mSolutions[ SolutionNewNotStim.mSolutions.size()-1 ];
       
        value1 = monodomain_pde.ComputeNonlinearSourceTermAtNode(node1, voltage);
        value2 = -(0 + monodomain_pde.GetIIonic(solutionSetNoStimT_05, voltage));

        TS_ASSERT_DELTA(value1, value2, 0.000001);
 
 #if 0
 
        // TS_ASSERT_DELTA(monodomain_pde.ComputeNonlinearSourceTermAtNode(node0, voltage),i_total,0.00001);

        // Check that we get the same result because ResetAsUnsolvedOdeSystem() has not yet been called.
        // TS_ASSERT_DELTA(monodomain_pde.ComputeNonlinearSourceTermAtNode(node0, voltage),i_total,0.00001);
        
        // Compare against known data with no stimulus at t = 0.5       
        i_stim  = 0.0;   
        v = -8.4507e+01;
        m = 1.6777e-03;
        h =  9.8328e-01; 
        j = 9.8950e-01;
        d = 2.9993e-03;  
        f = 1.0000e-00;
        x = 5.6001e-03;
        caI = 1.9926e-04;
        
        std::vector<double> solutionSetNoStimT_05;
        solutionSetNoStimT_05.push_back(h);
        solutionSetNoStimT_05.push_back(j);
        solutionSetNoStimT_05.push_back(m);
        solutionSetNoStimT_05.push_back(caI);
        solutionSetNoStimT_05.push_back(v);
        solutionSetNoStimT_05.push_back(d);
        solutionSetNoStimT_05.push_back(f);
        solutionSetNoStimT_05.push_back(x);
                              
        // remark: the 'voltage' used is the resting membrane potential -84.5mV                      
        i_ionic = monodomain_pde.GetIIonic(solutionSetNoStimT_05,voltage);
        i_total = i_stim + i_ionic;
         
        //TS_ASSERT_DELTA(monodomain_pde.ComputeNonlinearSourceTermAtNode(node1, voltage),i_total,0.00000001);
  
   
        // Reset       
        monodomain_pde.ResetAsUnsolvedOdeSystem();
                       
        voltage = -4.4120e+01; // reset voltage to value at t = 0.5

        //Comparing data against data for known solution: t = 1
        v = 4.2094e+01;
        m = 9.9982e-01;
        h = 9.4396e-02; 
        j = 8.7622e-01; 
        d = 2.1953e-02;
        f = 9.9889e-01;   
        x = 6.8886e-03; 
        caI = 1.9979e-04;
        
        std::vector<double> solutionSetStimT_1; 
        solutionSetStimT_1.push_back(h);
        solutionSetStimT_1.push_back(j);
        solutionSetStimT_1.push_back(m); 
        solutionSetStimT_1.push_back(caI);
        solutionSetStimT_1.push_back(v);
        solutionSetStimT_1.push_back(d);
        solutionSetStimT_1.push_back(f);
        solutionSetStimT_1.push_back(x);
        
        magnitudeOfStimulus = 0; //at time t=1, stimulus is zero
        i_stim = magnitudeOfStimulus;
                
        i_ionic = monodomain_pde.GetIIonic(solutionSetStimT_1, voltage);
        i_total = i_stim + i_ionic; 
         
     //   TS_ASSERT_DELTA(monodomain_pde.ComputeNonlinearSourceTermAtNode(node0, voltage),i_total,0.01);

        // Check that we get the same result because ResetAsUnsolvedOdeSystem() has not yet been called.
     //   TS_ASSERT_DELTA(monodomain_pde.ComputeNonlinearSourceTermAtNode(node0, voltage),i_total,0.01); 

        
        // Compare against known data with no stimulus at t=1
        voltage = -84.512;
        i_stim  = 0.0;   
        v = -8.4512e+01;
        m =  1.6761e-03;
        h = 9.8327e-01; 
        j =  9.8950e-01;
        d = 2.9986e-03;
        f = 1.0000e-00;
        x = 5.6003e-03;
        caI = 1.9854e-04;
                              
        std::vector<double> solutionSetNoStimT_1; 
        solutionSetNoStimT_1.push_back(h);
        solutionSetNoStimT_1.push_back(j);
        solutionSetNoStimT_1.push_back(m);
        solutionSetNoStimT_1.push_back(caI);
        solutionSetNoStimT_1.push_back(v);
        solutionSetNoStimT_1.push_back(d);
        solutionSetNoStimT_1.push_back(f);
        solutionSetNoStimT_1.push_back(x);
        
        i_ionic = monodomain_pde.GetIIonic(solutionSetNoStimT_1, voltage);
        i_total = i_stim + i_ionic;
        
    //    TS_ASSERT_DELTA(monodomain_pde.ComputeNonlinearSourceTermAtNode(node1, voltage),i_total,0.000000001);


#endif

    }
    
     
};




#endif //_TESTMONODOMAINPDE_HPP_
