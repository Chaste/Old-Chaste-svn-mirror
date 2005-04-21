#ifndef _TESTODESOLVERFORLR91_HPP_
#define _TESTODESOLVERFORLR91_HPP_


#include <cxxtest/TestSuite.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "LR91Model.hpp"
#include "AbstractStimulusFunction.hpp"
#include "InitialStimulus.hpp"


class TestOdeSolverForLR91 : public CxxTest::TestSuite
{
    public:
    
    // Test Ode Solver for LR91
    void testOdeSolverForLR91(void)
    {
    //Tests that setting and getting stimulus works -test passed
//        double time = 0.8;
//        double magnitudeOfStimulus = 80.0;        
//        AbstractStimulusFunction *pStimulus = new InitialStimulus(magnitudeOfStimulus);
//        double stimulusValue = pStimulus->GetStimulus(time);
//        std::cout<<"Stimulus value is "<< stimulusValue<< std::endl;
//        
//        // test 1
        double time = 0.2;
        double magnitudeOfStimulus = 80.0;        
        AbstractStimulusFunction *pStimulus = new InitialStimulus(magnitudeOfStimulus);
        LR91OdeFun lr91OdeFun(pStimulus);
        double voltage = -40.0;
        double m = 0.1;
        double h = 0.03;
        double j = 0.7;
        double d = 0.05;
        double f = 0.5;
        double x = 0.1;
        double caI = 0.001;
          std::vector<double> Y(8);
          Y[0] = voltage;
          Y[1] = m;
          Y[2] = h;
          Y[3] = j;
          Y[4] = d;
          Y[5] = f;
          Y[6] = x;
          Y[7] = caI;
        //std::vector<double> Result = lr91OdeFun.EvaluateYDerivatives(time, Y);
        
//        
//       //test 2 
//        voltage = -40.0;
//        m = 0.1;
//        h = 0.03;
//        j = 0.7;
//        d = 0.05;
//        f = 0.5;
//        x = 0.1;
//        caI = 0.001;
//        double magnitudeOfStimulus = 80.0;        
//        AbstractStimulusFunction *pStimulus = new InitialStimulus(magnitudeOfStimulus);
//        
//        LR91Model *pLR91Model;
//        pLR91Model = new LR91Model(v, m, double h, double j, double d, 
//                                    double f, double x, double caI, AbstractStimulusFunction stimulus stimulus);
//        std::vector<double> Result = pLR91Model->Solve();
//        output to matlab file
    }
    
};



#endif //_TESTODESOLVERFORLR91_HPP_
