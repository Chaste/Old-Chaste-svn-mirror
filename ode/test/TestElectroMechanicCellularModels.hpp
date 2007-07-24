#ifndef TESTELECTROMECHANICCALELLULARMODELS_HPP_
#define TESTELECTROMECHANICCALELLULARMODELS_HPP_

#include <cxxtest/TestSuite.h>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include "NHSCellularMechanicsOdeSystem.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "InitialStimulus.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"

double MyLam(double t, double endTime)
{
    double scaled_t = 10*t/endTime - 5;
    return 0.9 + 0.1*tanh(scaled_t);
}
double MyLamDeriv(double t, double endTime)
{
    double scaled_t = 10*t/endTime - 5;
    return (1.0/endTime)*(1-tanh(scaled_t)*tanh(scaled_t));
}


class TestElectroMechanicCellularModels : public CxxTest::TestSuite
{
public:
    void TestSingleCellWithSpecifiedLambda()
    {
        // setup               
        double magnitude = -25.5;
        double duration  =   2.0;  // ms
        double when      =   0.0;  // ms
        InitialStimulus stimulus(magnitude, duration, when);
        
        double end_time = 1.0; 
        double time_step = 0.001; 
               
        EulerIvpOdeSolver solver;
        LuoRudyIModel1991OdeSystem electrophys_model(&solver, time_step, &stimulus);
        NHSCellularMechanicsOdeSystem cellmech_model;
        
        // find out if electrophys model has CaTrop
        unsigned Ca_i_index = electrophys_model.GetStateVariableNumberByName("CaI");
        bool has_Ca_trop = false;
        unsigned Ca_trop_index;
        
        try
        {
            Ca_trop_index = electrophys_model.GetStateVariableNumberByName("CaTrop");
            has_Ca_trop = true;
        }
        catch(Exception& e)
        {
            has_Ca_trop = false;
        }

        std::vector<double> times;
        std::vector<double> active_tensions;
                
        // time loop
        for(double current_time = 0; current_time<end_time; current_time+=time_step)
        {
            // solve electrophys model
            electrophys_model.Compute(current_time, current_time+time_step);
            
            // get CaI 
            double Ca_I = electrophys_model.rGetStateVariables()[Ca_i_index];
            
            // set lam, dlam_dt and Ca_i in the cellular mechanics model
            double lam = MyLam(current_time, end_time);
            double dlam_dt = MyLamDeriv(current_time, end_time);
            cellmech_model.SetLambda1DerivativeAndCalciumI(lam, dlam_dt, Ca_I);
            
            // solve the cellular mechanics model
            solver.SolveAndUpdateStateVariable(&cellmech_model, current_time, current_time+time_step, time_step);
                    
            // IF electrophy model has CaTrop, get CaTrop from
            // cellular mechanics models and send to electrophy model
            if(has_Ca_trop)
            {
                std::vector<double>& electrophys_model_state_vars = electrophys_model.rGetStateVariables();
                electrophys_model_state_vars[Ca_trop_index] = cellmech_model.GetCalciumTroponinValue();
            }
            
            std::cout << current_time << " " << cellmech_model.GetActiveTension() << "\n";

//            times.push_back(current_time);
//            active_tensions.push_back( cellmech_model.GetActiveTension() );
        }
        
    }
};
#endif /*TESTELECTROMECHANICCALELLULARMODELS_HPP_*/
