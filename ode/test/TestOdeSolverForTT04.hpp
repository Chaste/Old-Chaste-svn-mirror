#ifndef _TESTODESOLVERFORTT04_HPP_
#define _TESTODESOLVERFORTT04_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "AbstractStimulusFunction.hpp"
#include "InitialStimulus.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "RungeKutta2IvpOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "AdamsBashforthIvpOdeSolver.hpp"
#include "OdeSolution.hpp"
#include "ColumnDataWriter.hpp"
#include "TenTusscherModel2004OdeSystem.hpp"

/**
 * \todo NOT WORKING AS YET - HAVE TO CHECK CellML model.....
 */

class TestOdeSolverForTT04 : public CxxTest::TestSuite
{
    public:
    
    // Test Ode Solver for TT04
    void TestOdeSolverForTT04WithNoStimulus(void)
    {
        /*
         * Set initial conditions and magnitude of stimulus
         * 
         */
        double calcium_dynamics_Ca_i = 0.0002;
        double calcium_dynamics_Ca_SR = 0.2;
        double calcium_dynamics_g = 1.0;
        double fast_sodium_current_h_gate_h = 0.75;
        double fast_sodium_current_j_gate_j = 0.75;
        double fast_sodium_current_m_gate_m = 0.0;
        double L_type_calcium_current_d_gate_d = 0.0;
        double L_type_calcium_current_f_Ca_gate_f_Ca = 1.0;
        double L_type_calcium_current_f_gate_f = 1.0;
        double membrane_V = -86.2;
        double potassium_dynamics_K_i = 138.3;
        double rapid_delayed_rectifier_current_X_r1_gate_X_r1 = 0.0;
        double rapid_delayed_rectifier_current_X_r2_gate_X_r2 = 1.0;
        double slow_delayed_rectifier_current_X_s_gate_X_s = 0.0;
        double sodium_dynamics_Na_i = 11.6;
        double transient_outward_current_r_gate_r = 0.0;
        double transient_outward_current_s_gate_s = 1.0;      
        
        double magnitudeOfStimulus = 0.0;  
        double durationOfStimulus  = 1.0 ;  // ms                     
        
        /*
         * Collect initial data in a vector
         * 
         */  
        std::vector<double> initialConditions;
        initialConditions.push_back(calcium_dynamics_Ca_i);
        initialConditions.push_back(calcium_dynamics_Ca_SR);
        initialConditions.push_back(calcium_dynamics_g);
        initialConditions.push_back(fast_sodium_current_h_gate_h);
        initialConditions.push_back(fast_sodium_current_j_gate_j);
        initialConditions.push_back(fast_sodium_current_m_gate_m);
        initialConditions.push_back(L_type_calcium_current_d_gate_d);
        initialConditions.push_back(L_type_calcium_current_f_Ca_gate_f_Ca);
        initialConditions.push_back(L_type_calcium_current_f_gate_f);
        initialConditions.push_back(membrane_V);
        initialConditions.push_back(potassium_dynamics_K_i);
        initialConditions.push_back(rapid_delayed_rectifier_current_X_r1_gate_X_r1);
        initialConditions.push_back(rapid_delayed_rectifier_current_X_r2_gate_X_r2);
        initialConditions.push_back(slow_delayed_rectifier_current_X_s_gate_X_s);
        initialConditions.push_back(sodium_dynamics_Na_i);
        initialConditions.push_back(transient_outward_current_r_gate_r);
        initialConditions.push_back(transient_outward_current_s_gate_s);
        /*
         * Choose function for stimulus
         */             
        InitialStimulus stimulus(magnitudeOfStimulus, durationOfStimulus); 
        
        /*
         * Instantiate the Luo-Rudy model: need to pass initial data and stimulus function
         */        
        TenTusscherModel2004OdeSystem tt04_ode_system(&stimulus);
        
        /*
         * Choose an ode solver
         */      
        EulerIvpOdeSolver solver;
        
        /*
         * Solve 
         */
        double startTime = 0.0;
        double endTime = 10.0;
        double timeStep = 0.01;             
                
        OdeSolution solution = solver.Solve(&tt04_ode_system, startTime, endTime, timeStep, initialConditions);
        
              
        /*
         * Write data to a file NewLR91.dat using ColumnDataWriter
         */                                                           
		        
        ColumnDataWriter *mpTestWriter;
        mpTestWriter = new ColumnDataWriter("testoutput","TT04Result");
        mpTestWriter->DefineFixedDimension("Time","ms", solution.mSolutions.size());
        int time_var_id = mpTestWriter->DefineVariable("Time","ms");
//        int a_var_id = mpTestWriter->DefineVariable("calcium_dynamics_Ca_i","microMol");/* Might be wrong units but that's what's in the paper */
//        int b_var_id = mpTestWriter->DefineVariable("calcium_dynamics_Ca_SR","mMol");
//        int c_var_id = mpTestWriter->DefineVariable("calcium_dynamics_g","mMol");
        int d_var_id = mpTestWriter->DefineVariable("fast_sodium_current_h_gate_h","milliamperes");
        int e_var_id = mpTestWriter->DefineVariable("fast_sodium_current_j_gate_j","milliamperes");
        int f_var_id = mpTestWriter->DefineVariable("fast_sodium_current_m_gate_m","milliamperes");
        int g_var_id = mpTestWriter->DefineVariable("L_type_calcium_current_d_gate_d","milliamperes");
        int h_var_id = mpTestWriter->DefineVariable("L_type_calcium_current_f_Ca_gate_f_Ca","milliamperes");
        int i_var_id = mpTestWriter->DefineVariable("L_type_calcium_current_f_gate_f","milliamperes");
//        int j_var_id = mpTestWriter->DefineVariable("membrane_V","millivolts");
//        int k_var_id = mpTestWriter->DefineVariable("potassium_dynamics_K_i","mMol");
        int l_var_id = mpTestWriter->DefineVariable("rapid_delayed_rectifier_current_X_r1_gate_X_r1","milliamperes");
        int m_var_id = mpTestWriter->DefineVariable("rapid_delayed_rectifier_current_X_r2_gate_X_r2","milliamperes");
        int n_var_id = mpTestWriter->DefineVariable("slow_delayed_rectifier_current_X_s_gate_X_s","milliamperes");
//        int o_var_id = mpTestWriter->DefineVariable("sodium_dynamics_Na_i","mMol");
        int p_var_id = mpTestWriter->DefineVariable("transient_outward_current_r_gate_r","milliamperes");
        int q_var_id = mpTestWriter->DefineVariable("transient_outward_current_s_gate_s","milliamperes");
        mpTestWriter->EndDefineMode();
				
        for (int i = 0; i < solution.mSolutions.size(); i++) 
        {
            mpTestWriter->PutVariable(time_var_id, solution.mTime[i], i);
//            mpTestWriter->PutVariable(a_var_id, solution.mSolutions[i][0], i);
//            mpTestWriter->PutVariable(b_var_id, solution.mSolutions[i][1], i);
//            mpTestWriter->PutVariable(c_var_id, solution.mSolutions[i][2], i);
            mpTestWriter->PutVariable(d_var_id, solution.mSolutions[i][3], i);
            mpTestWriter->PutVariable(e_var_id, solution.mSolutions[i][4], i);
            mpTestWriter->PutVariable(f_var_id, solution.mSolutions[i][5], i);
            mpTestWriter->PutVariable(g_var_id, solution.mSolutions[i][6], i);
            mpTestWriter->PutVariable(h_var_id, solution.mSolutions[i][7], i);     
            mpTestWriter->PutVariable(i_var_id, solution.mSolutions[i][8], i);
//            mpTestWriter->PutVariable(j_var_id, solution.mSolutions[i][9], i);
//            mpTestWriter->PutVariable(k_var_id, solution.mSolutions[i][10], i);
            mpTestWriter->PutVariable(l_var_id, solution.mSolutions[i][11], i);
            mpTestWriter->PutVariable(m_var_id, solution.mSolutions[i][12], i);
            mpTestWriter->PutVariable(n_var_id, solution.mSolutions[i][13], i);
//            mpTestWriter->PutVariable(o_var_id, solution.mSolutions[i][14], i);
            mpTestWriter->PutVariable(p_var_id, solution.mSolutions[i][15], i);  
            mpTestWriter->PutVariable(q_var_id, solution.mSolutions[i][16], i);       
        }
        
        delete mpTestWriter;
        
        
//        //read in good data file and compare line by line
//        std::ifstream testfile("testoutput/TT04Result.dat",std::ios::in);
//        std::ifstream goodfile("ode/test/data/TT04ResultGood.dat",std::ios::in);
//        std::string teststring;
//        std::string goodstring;
//        while(getline(testfile, teststring))
//        {
//              getline(goodfile,goodstring);
//              TS_ASSERT_EQUALS(teststring,goodstring);
//        }
//        testfile.close();
//        goodfile.close();

    }	

};



#endif //_TESTODESOLVERFORLR91_HPP_
