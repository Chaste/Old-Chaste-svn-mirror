#ifndef TESTHODGKINHUXLEYSQUIDAXON1952ORIGINALODESYSTEM_HPP_
#define TESTHODGKINHUXLEYSQUIDAXON1952ORIGINALODESYSTEM_HPP_

#include <cxxtest/TestSuite.h>
//#include <iostream>
#include <vector>

#include "InitialStimulus.hpp"
#include "HodgkinHuxleySquidAxon1952OriginalOdeSystem.hpp"

class TestHodgkinHuxleySquidAxon1952OriginalOdeSystem: public CxxTest::TestSuite
{
    public:
   
    void testHHModelAtSingularities()
    {
         /*
         * Set stimulus
         */   
        double magnitude_stimulus = 0.0;  // uA/cm2
        double duration_stimulus = 0.;  // ms
        double start_stimulus = 0.0;   // ms
        InitialStimulus stimulus(magnitude_stimulus,
                                 duration_stimulus,
                                 start_stimulus);

        HodgkinHuxleySquidAxon1952OriginalOdeSystem hh52_ode_system(NULL, &stimulus, 0.01);
        
        double v_singularity[2];
        v_singularity[0]=-65;
        v_singularity[1]=-50;
        
        for (int i=0; i<2; i++)
        {
            
            std::vector<double> yleft;
           
            //mVariableNames.push_back("V");
            //mVariableUnits.push_back("mV");
            yleft.push_back(v_singularity[i]+0.1);
            
            //mVariableNames.push_back("n");
            //mVariableUnits.push_back("");
            yleft.push_back(0.325);
            
            //mVariableNames.push_back("h");
            //mVariableUnits.push_back("");
            yleft.push_back(0.6);
            
            //mVariableNames.push_back("m");
            //mVariableUnits.push_back("");
            yleft.push_back(0.05);
                 
        
            std::vector<double> rhsleft =  hh52_ode_system.EvaluateYDerivatives (0.0, yleft);
        
            std::vector<double> yright;
           
            //mVariableNames.push_back("V");
            //mVariableUnits.push_back("mV");
            yright.push_back(v_singularity[i]-0.1);
            
            //mVariableNames.push_back("n");
            //mVariableUnits.push_back("");
            yright.push_back(0.325);
            
            //mVariableNames.push_back("h");
            //mVariableUnits.push_back("");
            yright.push_back(0.6);
            
            //mVariableNames.push_back("m");
            //mVariableUnits.push_back("");
            yright.push_back(0.05);
                 
            std::vector<double> rhsright = hh52_ode_system.EvaluateYDerivatives (0.0, yright);
        
        
            std::vector<double> y_at_singularity;
           
            //mVariableNames.push_back("V");
            //mVariableUnits.push_back("mV");
            y_at_singularity.push_back(v_singularity[i]);
            
            //mVariableNames.push_back("n");
            //mVariableUnits.push_back("");
            y_at_singularity.push_back(0.325);
            
            //mVariableNames.push_back("h");
            //mVariableUnits.push_back("");
            y_at_singularity.push_back(0.6);
            
            //mVariableNames.push_back("m");
            //mVariableUnits.push_back("");
            y_at_singularity.push_back(0.05);
                 
            std::vector<double> rhs_at_singularity = hh52_ode_system.EvaluateYDerivatives (0.0, y_at_singularity);
            
            for (int j=0; j<4; j++)
            {
                // std::cout << j << std::endl;
                TS_ASSERT_DELTA((rhsright[j]+rhsleft[j])/2, rhs_at_singularity[j], 0.1);
            }
        }
    }
};

#endif /*TESTHODGKINHUXLEYSQUIDAXON1952ORIGINALODESYSTEM_HPP_*/
