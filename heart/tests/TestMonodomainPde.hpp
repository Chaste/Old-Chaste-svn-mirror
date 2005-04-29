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


#include <cxxtest/TestSuite.h>


class TestMonodomainPde : public CxxTest::TestSuite
{
    public:
    
    void testMonodomainPdeConstructor( void )
    {
        int num_nodes=2;
        
        Node<1> node0(0,true,0);
        Node<1> node1(1,true,0);
        
        double big_time_step=1;
        double small_time_step=0.01;
        AbstractIvpOdeSolver *pMySolver = new EulerIvpOdeSolver();
        MonodomainPde monodomain_pde(num_nodes, big_time_step, pMySolver, small_time_step);
        
        double voltage = -9999; // This voltage will be ignored
        double m = 0.0017;
        double h = 0.9833;
        double j = 0.9895;
        double d = 0.003;
        double f = 1;
        double x = 0.0056;
        double caI = 0.0002;
        double magnitudeOfStimulus = -80.0;  
        double durationOfStimulus  = 0.5 ;
                  
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
        //int node_index = 0;     
        monodomain_pde.SetStimulusFunctionAtNode(0, pStimulus);
        
        voltage = -84.5;
   

        double i_stim = 0.0; // Set i_stim to zero by default
        
        double i_ionic;
        double i_total;        
                  
        // __________________________________________________________________
        /*
         * If want values for t=1 ms solving odes with stimulus, use data below
         */
        double v = 4.2094e+01;
        m = 9.9982e-01;
        h = 9.4396e-02; 
        j = 8.7622e-01; 
        d = 2.1953e-02;
        f = 9.9889e-01;
        x = 6.8886e-03;
        caI = 1.9979e-04;
//        
        // these are values after 0.01 ms
//        double v = -8.3700e+01;
//       m = 1.6645e-03;
//        h =  9.8330e-01 ; 
//        j =9.8950e-01 ;
//        d = 3.0000e-03;
//        f =  1.0000e-00;
//        x = 5.6000e-03 ;
//        caI = 1.9998e-04;
        
        std::vector<double> solutionSet;
        solutionSet.push_back(h);
        solutionSet.push_back(j);
        solutionSet.push_back(m);
        solutionSet.push_back(caI);
        solutionSet.push_back(v);
        solutionSet.push_back(d);
        solutionSet.push_back(f);
        solutionSet.push_back(x);
        
        magnitudeOfStimulus = 0; //at time t=1, stimulus is zero
        i_stim = magnitudeOfStimulus;
        
        i_ionic = monodomain_pde.GetIIonic(solutionSet, voltage);
        i_total = i_stim + i_ionic;
        
        TS_ASSERT_DELTA(monodomain_pde.ComputeNonlinearSourceTermAtNode(node0, voltage),i_total,0.001);

        // Check that we get the same result because NextStep was not called.
        TS_ASSERT_DELTA(monodomain_pde.ComputeNonlinearSourceTermAtNode(node0, voltage),i_total,0.001);
        
         /*
         * If want values for t=1 ms with NO stimulus, use data below
         */       
         i_stim  = 0.0;   
         v = -8.4512e+01;
         m =  1.6761e-03;
         h = 9.8327e-01; 
         j =  9.8950e-01;
         d = 2.9986e-03;
         f = 1.0000e-00;
         x = 5.6003e-03;
         caI = 1.9854e-04;
                              
         std::vector<double> solutionSet2;
         solutionSet2.push_back(h);
         solutionSet2.push_back(j);
         solutionSet2.push_back(m);
         solutionSet2.push_back(caI);
         solutionSet2.push_back(v);
         solutionSet2.push_back(d);
         solutionSet2.push_back(f);
         solutionSet2.push_back(x);
                               
         i_ionic = monodomain_pde.GetIIonic(solutionSet2,voltage);
         i_total = i_stim + i_ionic;
        
         TS_ASSERT_DELTA(monodomain_pde.ComputeNonlinearSourceTermAtNode(node1, voltage),i_total,0.001);
     }
     
     
     //void testReset
     
     //-> in new test file TestDg0MonodomainAssembler
     //   test setup
     //   test update g first
     //   test update v first
     //   test with init stim and with reg stim with >1 node
};




#endif //_TESTMONODOMAINPDE_HPP_
