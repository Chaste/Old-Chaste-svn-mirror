/*

Copyright (C) University of Oxford, 2005-2009

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
#include <cmath>

#include "BackwardEulerLuoRudyIModel1991.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "PdeConvergenceTester.hpp"
#include "SpaceConvergenceTester.hpp"
#include "KspConvergenceTester.hpp"
#include "OdeConvergenceTester.hpp"
#include "OdePdeConvergenceTester.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "BackwardEulerNobleVargheseKohlNoble1998.hpp"
#include "NobleVargheseKohlNoble1998.hpp"
#include "NobleVargheseKohlNoble1998Optimised.hpp"


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
        HeartConfig::Instance()->SetUseAbsoluteTolerance(5e-4);
        //pTester->SetKspAbsoluteTolerance(5e-4);

        pTester->Converge("Automated_test");
        TS_ASSERT(pTester->Converged);
    }

    void ConvergeInVarious(StimulusType stimulusType)
    {
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetKSPPreconditioner(), "bjacobi");
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetKSPSolver(), "cg");
        HeartConfig::Instance()->SetKSPPreconditioner("jacobi");
        HeartConfig::Instance()->SetKSPSolver("gmres");
        {
            std::cout << "PdeConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1, 2>\n";
            PdeConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1, 2> tester;
            RunConvergenceTester(&tester, stimulusType);
            TS_ASSERT_DELTA(tester.PdeTimeStep, 5.0e-3, 1e-10);
        }

        {
            std::cout << "SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1, 2>\n";
            //Block Jacobi with CG can detect zero pivots in a 1-D convergence test
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
            TS_ASSERT_EQUALS(HeartConfig::Instance()->GetAbsoluteTolerance(), 5.0e-4);
            KspConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1, 2> tester;
            RunConvergenceTester(&tester, stimulusType);
            if (stimulusType != REGION)
            {
                TS_ASSERT_DELTA(HeartConfig::Instance()->GetAbsoluteTolerance(), 1e-3, 1e-10);
            }
            else
            {
                TS_ASSERT_DELTA(HeartConfig::Instance()->GetAbsoluteTolerance(), 1e-4, 1e-10);
            }
            //See above - we've fiddled with HeartConfig...
            HeartConfig::Instance()->SetUseAbsoluteTolerance(5.0e-4);
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
        //Put the KSP defaults back (see above)
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetKSPPreconditioner(), "jacobi");
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetKSPSolver(), "gmres");
        HeartConfig::Instance()->SetKSPPreconditioner("bjacobi");
        HeartConfig::Instance()->SetKSPSolver("cg");
    }

public:

    void joeTestFullActionPotential() throw(Exception)
    {
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1, 2> tester;
        tester.SimulateFullActionPotential=true;
        //Time steps are okay for giving a sensible upstroke
        tester.PdeTimeStep=0.1;
        tester.OdeTimeStep=0.1;
        
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.IsConverged());
        
        ///Note that long plateau phase will force convergence to happen earlier 
        TS_ASSERT_EQUALS(tester.MeshNum, 4u);
        
        TS_ASSERT_DELTA(329.0, tester.Apd90FirstQn, 1.5);
        TS_ASSERT_DELTA(329.0, tester.Apd90ThirdQn, 1.5);
        TS_ASSERT_DELTA(0.0588, tester.ConductionVelocity, 1e-3);
    }

    void TestStimulatePlanein1D() throw(Exception)
    {
        ConvergeInVarious(PLANE);
    }

    void TestStimulateRegionin1D() throw(Exception)
    {
        ConvergeInVarious(REGION);
    }
    

    //Current test takes about 20 mins.
    //This is much longer (1 hour?) with default ksp
    void Test2DSpaceSymmLq() throw(Exception)
    {
        HeartConfig::Instance()->SetKSPSolver("symmlq");
        HeartConfig::Instance()->SetKSPPreconditioner("bjacobi");
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<2>, 2, 2> tester;
        //tester.SetKspAbsoluteTolerance(1e-3);
         HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-3);
        
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 5u);
        TS_ASSERT(tester.IsConverged());
        //TS_ASSERT_EQUALS(tester.GetMeshNum(), 4);
        //TS_ASSERT_DELTA(tester.GetSpaceStep(), 0.0023 /*cm*/, 1e-4 /*Allowed error*/);
        HeartConfig::Instance()->Reset();
    }

    void Test2DSpaceWithRegionStimulus() throw(Exception)
    {
        HeartConfig::Instance()->SetKSPSolver("symmlq");
        HeartConfig::Instance()->SetKSPPreconditioner("bjacobi");
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<2>, 2, 2> tester;
        tester.Stimulus = REGION;
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-3);
        
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 6u);
        HeartConfig::Instance()->Reset();
    }



    //Currently takes about 3 minutes to do mesh0 and mesh1
    void Test3DSpace() throw(Exception)
    {
        HeartConfig::Instance()->SetKSPSolver("symmlq");
        HeartConfig::Instance()->SetKSPPreconditioner("bjacobi");
        SpaceConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<3>, 3, 2> tester;
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-3);
        
        tester.RelativeConvergenceCriterion=4e-2;//Just to prove the thing works
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 1u); ///Just to prove the thing works
        HeartConfig::Instance()->Reset();
    }

    void TestSpaceConvergencein1DWithBackwardN98() throw(Exception)
    {
        SpaceConvergenceTester<BackwardEulerNobleVargheseKohlNoble1998,  MonodomainProblem<1>, 1, 1> tester;
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_EQUALS(tester.MeshNum, 5u);
        TS_ASSERT_LESS_THAN(tester.LastDifference, 1.68417e-05);
    }

    void TestOdeConvergencein1DWithBackwardN98() throw(Exception)
    {
        OdeConvergenceTester<BackwardEulerNobleVargheseKohlNoble1998,  MonodomainProblem<1>, 1, 1> tester;
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_DELTA(tester.OdeTimeStep, 0.005, 1e-10);
    }

    void TestOdePdeConvergencein1DWithBackwardN98() throw(Exception)
    {
        OdePdeConvergenceTester<BackwardEulerNobleVargheseKohlNoble1998,  MonodomainProblem<1>, 1, 1> tester;
        tester.NeumannStimulus = 5000;
        tester.Stimulus = NEUMANN;
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_DELTA(tester.OdeTimeStep, 0.005, 1e-10);
        TS_ASSERT_DELTA(tester.PdeTimeStep, 0.005, 1e-10);
    }

    void TestOdePdeConvergencein1DWithForwardLookupN98() throw(Exception)
    {
        OdePdeConvergenceTester<CML_noble_varghese_kohl_noble_1998_basic_pe_lut,  MonodomainProblem<1>, 1, 1> tester;
        tester.NeumannStimulus = 5000;
        tester.Stimulus = NEUMANN;
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_DELTA(tester.OdeTimeStep, 0.0025, 1e-10);
        TS_ASSERT_DELTA(tester.PdeTimeStep, 0.0025, 1e-10);
    }
    void TestOdePdeConvergencein1DWithForwardBasicN98() throw(Exception)
    {
        OdePdeConvergenceTester<CML_noble_varghese_kohl_noble_1998_basic,  MonodomainProblem<1>, 1, 1> tester;
        tester.NeumannStimulus = 5000;
        tester.Stimulus = NEUMANN;
        tester.Converge(__FUNCTION__);
        TS_ASSERT(tester.Converged);
        TS_ASSERT_DELTA(tester.OdeTimeStep, 0.0025, 1e-10);
        TS_ASSERT_DELTA(tester.PdeTimeStep, 0.0025, 1e-10);
    }
};

#endif /*TESTCONVERGENCENIGHTLY_HPP_*/
