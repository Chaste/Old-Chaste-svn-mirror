/*

Copyright (C) University of Oxford, 2008

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


#ifndef TESTCONVERGENCENIGHTLY_HPP_
#define TESTCONVERGENCENIGHTLY_HPP_

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
#include "PetscSetupAndFinalize.hpp"


class TestConvergenceNightly : public CxxTest::TestSuite
{

public:

    void RunConvergenceTester(AbstractUntemplatedConvergenceTester *pTester, StimulusType stimulusType)
    {
        pTester->Stimulus = stimulusType;
        if ( stimulusType == REGION )
        {
            pTester->MeshNum = 6u;
        }

        pTester->Converge("Automated_test");
        TS_ASSERT(pTester->Converged);
    }

    void ConvergeInVarious(StimulusType stimulusType)
    {
        {
            std::cout << "PdeConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1, 2>\n";
            PdeConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1, 2> tester;
            RunConvergenceTester(&tester, stimulusType);
            TS_ASSERT_DELTA(tester.PdeTimeStep, 5.0e-3, 1e-10);
        }

        {
            std::cout << "SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1, 2>\n";
            SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1, 2> tester;
            RunConvergenceTester(&tester, stimulusType);
            if (stimulusType != REGION)
            {
                TS_ASSERT_EQUALS(tester.MeshNum, 5u);
            }
            else
            {
                TS_ASSERT_EQUALS(tester.MeshNum, 6u);
            }
        }

        {
            std::cout << "KspConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1, 2>\n";
            KspConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1, 2> tester;
            RunConvergenceTester(&tester, stimulusType);
            if (stimulusType != REGION)
            {
                TS_ASSERT_DELTA(tester.GetKspRelativeTolerance(), 1e-5, 1e-10);
            }
            else
            {
                TS_ASSERT_DELTA(tester.GetKspRelativeTolerance(), 1e-6, 1e-10);
            }
        }

        {
            std::cout << "OdeConvergenceTester<LuoRudyIModel1991OdeSystem, BidomainProblem<1>, 1, 2>\n";
            OdeConvergenceTester<LuoRudyIModel1991OdeSystem, BidomainProblem<1>, 1, 2> tester;
            RunConvergenceTester(&tester, stimulusType);
            TS_ASSERT_DELTA(tester.OdeTimeStep, 0.0025, 1e-10);
        }

        {
            std::cout << "OdeConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1, 2>\n";
            OdeConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1, 2> tester;
            tester.PdeTimeStep=0.01;
            RunConvergenceTester(&tester, stimulusType);
            TS_ASSERT_DELTA(tester.OdeTimeStep, 0.0025, 1e-10);
        }

    }

public:

    void TestStimulatePlanein1D() throw(Exception)
    {
        ConvergeInVarious(PLANE);
    }

    void TestStimulateRegionin1D() throw(Exception)
    {
        ConvergeInVarious(REGION);
    }

    void TestVariousWithNeumannStimulus() throw(Exception)
    {
        ConvergeInVarious(NEUMANN);
    }


    //Current test takes about 20 mins.
    //This is much longer (1 hour?) with default ksp
    void Test2DSpaceSymmLq() throw(Exception)
    {
        HeartConfig::Instance()->SetKSPSolver("symmlq");
        HeartConfig::Instance()->SetKSPPreconditioner("bjacobi");
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<2>, 2, 2> tester;
        tester.SetKspRelativeTolerance(5e-8);
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 5u);
        TS_ASSERT(tester.IsConverged());
        //TS_ASSERT_EQUALS(tester.GetMeshNum(), 4);
        //TS_ASSERT_DELTA(tester.GetSpaceStep(), 0.0023 /*cm*/, 1e-4 /*Allowed error*/);
    }

    void Test2DSpaceSymmLqWithNeumannStimulus() throw(Exception)
    {
        HeartConfig::Instance()->SetKSPSolver("symmlq");
        HeartConfig::Instance()->SetKSPPreconditioner("bjacobi");
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<2>, 2, 2> tester;
        tester.SetKspRelativeTolerance(1e-9);
        tester.Stimulus = NEUMANN;
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 5u);
        TS_ASSERT(tester.IsConverged());
    }

    void Test2DSpaceWithRegionStimulus() throw(Exception)
    {
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<2>, 2, 2> tester;
        tester.Stimulus = REGION;
        tester.SetKspRelativeTolerance(1e-8);
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 6u);
    }

    void Test2DSpaceWithNeumannStimulus() throw(Exception)
    {
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<2>, 2, 2> tester;
        tester.Stimulus = NEUMANN;
        tester.SetKspRelativeTolerance(1e-8);
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 5u);
    }

    //Currently takes about 3 minutes to do mesh0 and mesh1
    void Test3DSpace() throw(Exception)
    {
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<3>, 3, 2> tester;
        tester.SetKspRelativeTolerance(1e-8);
        tester.RelativeConvergenceCriterion=4e-2;//Just to prove the thing works
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 1u); ///Just to prove the thing works
    }

    void Test3DSpaceWithNeumannStimulus() throw(Exception)
    {
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<3>, 3, 2> tester;
        tester.SetKspRelativeTolerance(1e-8);
        tester.RelativeConvergenceCriterion = 4e-2;//Just to prove the thing works
        tester.Stimulus = NEUMANN;
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 2u); ///Just to prove the thing works
    }
};

#endif /*TESTCONVERGENCENIGHTLY_HPP_*/
