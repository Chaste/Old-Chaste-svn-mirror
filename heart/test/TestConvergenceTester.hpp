/*
Copyright (C) University of Oxford, 2008

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with Chaste.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef TESTCONVERGENCETESTER_HPP_
#define TESTCONVERGENCETESTER_HPP_

#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "MonodomainProblem.hpp"
#include <petscvec.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>

#include "BackwardEulerLuoRudyIModel1991.hpp"
#include "PdeConvergenceTester.hpp"
#include "SpaceConvergenceTester.hpp"
#include "KspConvergenceTester.hpp"
#include "OdeConvergenceTester.hpp"

class TestConvergenceTester : public CxxTest::TestSuite
{   
public:
    void Test1DOdeTime() throw(Exception)
    {
        OdeConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1, 2> tester;
        tester.MeshNum=1;
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT(tester.IsConverged());
        TS_ASSERT_EQUALS(tester.OdeTimeStep, 0.0025); 
    }
    
    void Test1DPdeTime() throw(Exception)
    {
        PdeConvergenceTester<BackwardEulerLuoRudyIModel1991, MonodomainProblem<1>, 1, 1> tester;
        tester.MeshNum=1;
        tester.RelativeConvergenceCriterion=7e-4;
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.PdeTimeStep, 0.01); 
    }
    
    void Test1DPdeTimeRegion() throw(Exception)
    {
        PdeConvergenceTester<BackwardEulerLuoRudyIModel1991, MonodomainProblem<1>, 1, 1> tester;
        tester.MeshNum=1;
        tester.Stimulus=REGION;
        tester.RelativeConvergenceCriterion=5e-4;
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.PdeTimeStep, 0.01); 
    }
    
    void Test1DPdeTimeNeumann() throw(Exception)
    {
        PdeConvergenceTester<BackwardEulerLuoRudyIModel1991, MonodomainProblem<1>, 1, 1> tester;
        tester.MeshNum=1;
        tester.Stimulus=NEUMANN;
        tester.RelativeConvergenceCriterion=5e-4;
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.PdeTimeStep, 0.02); 
    }

    void Test1DSpace() throw(Exception)
    {
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, MonodomainProblem<1>, 1, 1> tester;
        tester.RelativeConvergenceCriterion=2e-2;
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.IsConverged());
        TS_ASSERT_EQUALS(tester.GetMeshNum(), 2); 
        TS_ASSERT_DELTA(tester.GetSpaceStep(), 0.0125, 1e-8); 
    }
    
    void Test2DOdeTime() throw(Exception)
    {
        OdeConvergenceTester<BackwardEulerLuoRudyIModel1991, MonodomainProblem<2>, 2, 1> tester;
        tester.MeshNum=0;
        tester.RelativeConvergenceCriterion=8e-5;
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.OdeTimeStep, 0.0025); 
    }
    
    void TestSpaceConvergencein1DWithRelativeTolerance() throw(Exception)
    {
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1, 2> tester;
        tester.SetKspRelativeTolerance(1e-4);
        TS_ASSERT_DELTA(tester.GetKspRelativeTolerance(), 1e-4, 1e-10);
        TS_ASSERT_THROWS_ANYTHING(tester.GetKspAbsoluteTolerance());
        tester.RelativeConvergenceCriterion=2e-2;
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.IsConverged());
        TS_ASSERT_EQUALS(tester.GetMeshNum(), 2); 
        TS_ASSERT_DELTA(tester.GetSpaceStep(), 0.0125, 1e-8);
        TS_ASSERT_LESS_THAN(tester.LastDifference, 3.1e-3);
    }
    
    void TestSpaceConvergencein1DWithAbsoluteTolerance() throw(Exception)
    {
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1, 2> tester;
        tester.SetKspAbsoluteTolerance(1e-4);
        TS_ASSERT_DELTA(tester.GetKspAbsoluteTolerance(), 1e-4, 1e-10);
        TS_ASSERT_THROWS_ANYTHING(tester.GetKspRelativeTolerance());
        tester.RelativeConvergenceCriterion=2e-2;
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.IsConverged());
        TS_ASSERT_EQUALS(tester.GetMeshNum(), 2); 
        TS_ASSERT_DELTA(tester.GetSpaceStep(), 0.0125, 1e-8);
        TS_ASSERT_LESS_THAN(tester.LastDifference, 2.9e-3);
     }
};

#endif /*TESTCONVERGENCETESTER_HPP_*/
