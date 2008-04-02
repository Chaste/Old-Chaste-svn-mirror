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

#ifndef TESTCONVERGENCENEUMANNSTIMULUS_HPP_
#define TESTCONVERGENCENEUMANNSTIMULUS_HPP_

#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "MonodomainProblem.hpp"
#include <petscvec.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>

#include "BackwardEulerLuoRudyIModel1991.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "PdeConvergenceTester.hpp"
#include "SpaceConvergenceTester.hpp"
#include "KspConvergenceTester.hpp"
#include "OdeConvergenceTester.hpp"
//#include "StimulusConvergenceTester.hpp"

class TestConvergenceNeumannStimulus : public CxxTest::TestSuite
{   
public:

    void TestConvergenceMonodomain1d()
    {
        PdeConvergenceTester<BackwardEulerLuoRudyIModel1991, MonodomainProblem<1>, 1, 1> tester;
        tester.Stimulus=NEUMANN;
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_DELTA(tester.PdeTimeStep, 5.0e-3, 1e-10);
    }

    void TestConvergenceBidomain1d()
    {
        PdeConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1, 2> tester;
        tester.Stimulus=NEUMANN;
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_DELTA(tester.PdeTimeStep, 5.0e-3, 1e-10);
    }

    void TestSpaceConvergence1d()
    {
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, MonodomainProblem<1>, 1, 1> tester;
        tester.Stimulus=NEUMANN;
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 5u); 
    }
    
    void TestSpaceConvergence2d()
    {
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, MonodomainProblem<2>, 2, 1> tester;
        tester.Stimulus=NEUMANN;
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 5u); 
    }
    
    void TestSpaceConvergence2dBidomain()
    {
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<2>, 2, 2> tester;
        tester.Stimulus=NEUMANN;
        tester.Converge();
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 5u); 
    }

};

#endif /*TESTCONVERGENCENEUMANNSTIMULUS_HPP_*/
