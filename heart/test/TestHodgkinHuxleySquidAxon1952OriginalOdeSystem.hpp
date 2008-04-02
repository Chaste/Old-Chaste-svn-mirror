/*
Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>.
*/

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

    void TestHHModelAtSingularities()
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
                                 
        HodgkinHuxleySquidAxon1952OriginalOdeSystem hh52_ode_system(NULL, 0.01, &stimulus);
        
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
            
            
            std::vector<double> rhsleft(yleft.size());
            hh52_ode_system.EvaluateYDerivatives (0.0, yleft, rhsleft);
            
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
            
            std::vector<double> rhsright(yright.size());
            hh52_ode_system.EvaluateYDerivatives (0.0, yright, rhsright);
            
            
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
            
            std::vector<double> rhs_at_singularity(y_at_singularity.size());
            hh52_ode_system.EvaluateYDerivatives (0.0, y_at_singularity, rhs_at_singularity);
            
            for (int j=0; j<4; j++)
            {
                // std::cout << j << std::endl;
                TS_ASSERT_DELTA((rhsright[j]+rhsleft[j])/2, rhs_at_singularity[j], 0.1);
            }
        }
    }
};

#endif /*TESTHODGKINHUXLEYSQUIDAXON1952ORIGINALODESYSTEM_HPP_*/
