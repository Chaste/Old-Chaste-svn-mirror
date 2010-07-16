/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


#ifndef _TESTIONICMODELSFAILING_HPP_
#define _TESTIONICMODELSFAILING_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <ctime>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "Exception.hpp"
#include "HeartConfig.hpp"

#include "ZeroStimulus.hpp"

#include "EulerIvpOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"

#include "TenTusscher2006Epi.hpp"
#include "TenTusscher2006EpiBackwardEuler.hpp"
#include "ArchiveLocationInfo.hpp"



class TestIonicModelsFailing : public CxxTest::TestSuite
{
public:


    void TestBackwardEulerBug() throw (Exception)
    {

        //These data come from Miguel's test (see #1339)        
        double dodgy_state_vars_array[19] = {-6.15475,0.00808679,0.284434,0.00633525,0.994096,0.0321343,0.402544,0.730188,0.856068,0.959682,0.998295,0.912007,0.02408,0.000115671,3.63196,0.00175538,0.937932,8.62141,136.891};
        std::vector<double> dodgy_state_vars(dodgy_state_vars_array,dodgy_state_vars_array+19);
        
        double step=0.1; // This was supposed to fail with step=0.05, but doesn't -- due to the accuracy of the initial state variables?
        
       // step=0.05; //Will pass REMOVE TO WORK ON #1339
        
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(step, step, step);
            
        boost::shared_ptr<ZeroStimulus> p_stimulus(new ZeroStimulus());
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        CellTenTusscher2006EpiFromCellMLBackwardEuler tt06_backward_euler(p_solver, p_stimulus);
        
        tt06_backward_euler.rGetStateVariables() = dodgy_state_vars;
        tt06_backward_euler.ComputeExceptVoltage(0.0, 3*step);
/*
The test from Miguel (might converge with a bit more work)
DEBUG: counter = 1, norm = 0.0076568
DEBUG: counter = 2, norm = 0.00494596
DEBUG: counter = 3, norm = 0.00198722
DEBUG: counter = 4, norm = 0.00556094
DEBUG: counter = 5, norm = 0.00538928
DEBUG: counter = 6, norm = 0.00952744
DEBUG: counter = 7, norm = 0.0232987
DEBUG: counter = 8, norm = 0.0144165
DEBUG: counter = 9, norm = 0.00606122
DEBUG: counter = 10, norm = 0.00353937
DEBUG: counter = 11, norm = 0.00196135
DEBUG: counter = 12, norm = 0.0148738
DEBUG: counter = 13, norm = 0.0230241
DEBUG: counter = 14, norm = 0.0108333
DEBUG: counter = 15, norm = 0.000441591
DEBUG: counter = 16, norm = 1.92038e-06
DEBUG: counter = 17, norm = 3.78497e-11
*/
    }
};


#endif //_TESTIONICMODELSFAILING_HPP_
