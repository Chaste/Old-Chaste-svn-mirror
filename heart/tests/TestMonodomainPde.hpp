#ifndef _TESTMONODOMAINPDE_HPP_
#define _TESTMONODOMAINPDE_HPP_

#include <cxxtest/TestSuite.h>

class TestMonodomainPde : public CxxTest::TestSuite
{
    public:
    
//    void testMonodomainPdeConstructor( void )
//    {
//        int num_nodes=2;
//        double big_time_step=1.0;
//        MonodomainPde monodomain_pde(num_nodes, big_time_step);
//        
//        double voltage = -9999; // This voltage will be ignored
//        double m = 0.0017;
//        double h = 0.9833;
//        double j = 0.9895;
//        double d = 0.003;
//        double f = 1;
//        double x = 0.0056;
//        double caI = 0.0002;
//        double magnitudeOfStimulus = -80.0;  
//        double durationOfStimulus  = 0.5 ;
//                  
//        std::vector<double> initialConditions;
//        initialConditions.push_back(h);
//        initialConditions.push_back(j);
//        initialConditions.push_back(m);
//        initialConditions.push_back(caI);
//        initialConditions.push_back(voltage);
//        initialConditions.push_back(d);
//        initialConditions.push_back(f);
//        initialConditions.push_back(x);
//        
//        MonodomainPde.SetUniversalInitialConditions(initialConditions);
//        
//        AbstractStimulusFunction *pStimulus = new InitialStimulus(magnitudeOfStimulus, durationOfStimulus);
//        int node_num = 0;     
//        MonodomainPde.SetStimulusFunctionAtNode(node_num, pStimulus);
//        
//        voltage = -84.5;
//        
//        mAtT = 9.9982e-01;
//        jAtT = 9.4396e-02;
//        
//         8.7622e-01  2.1953e-02  9.9889e-01  6.8886e-03  1.9979e-04
//         
//        plateau_potassium_current_i_Kp = plateau_potassium_current_g_Kp*plateau_potassium_current_Kp*(membrane_V-plateau_potassium_current_E_Kp);
//        
//        double i_ionic = pow( 
//        
//        
//        TS_ASSERT_DELTA(MonodomainPde.ComputeLinearSourceTermAtNode(node_num, voltage),
//                        i_ionic,
//                        0.000001);

};

#endif //_TESTMONODOMAINPDE_HPP_
