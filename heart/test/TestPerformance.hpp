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

#ifndef TESTPERFORMANCE_HPP_
#define TESTPERFORMANCE_HPP_

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
#include "PerformanceTester.hpp"

class TestPerformance : public CxxTest::TestSuite
{   
public:
    void TestPerf() throw(Exception)
    {
        PetscOptionsSetValue("-ksp_type", "symmlq");
        PetscOptionsSetValue("-pc_type", "bjacobi");
        PetscOptionsSetValue("-options_table", "");
        // write headings
        PerformanceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<2>, 2>::DisplayHeadings();
        EventHandler::Headings();

        // base line test
        PerformanceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<2>, 2> tester;
        tester.MeshNum=2;
        tester.SimTime=4.0;
        tester.Run();
        EventHandler::Report();

        // vary simulation time
        tester.SimTime=0.0025;
        tester.Run();
        EventHandler::Report();

        tester.SimTime=8.0;
        tester.Run();
        EventHandler::Report();
        
        tester.SimTime=4.0;
        
        // vary pde time step
        tester.PdeTimeStep=0.005;
        tester.Run();
        EventHandler::Report();
        
        tester.PdeTimeStep=0.01;
        tester.Run();
        EventHandler::Report();
                
        tester.PdeTimeStep=0.0025;        
        // vary ode time step
        tester.OdeTimeStep=0.0025/2;
        tester.Run();
        EventHandler::Report();       

        tester.OdeTimeStep=0.0025/4;
        tester.Run();
        EventHandler::Report();

        tester.OdeTimeStep=0.0025;
        
        // vary printing time step
        tester.PrintingTimeStep = 0.02;
        tester.Run();
        EventHandler::Report();
        
        tester.PrintingTimeStep = 0.01;
        tester.Run();
        EventHandler::Report();

        tester.PrintingTimeStep = 0.04;
        
        // vary mesh size
        tester.MeshNum++;
        tester.Run();
        EventHandler::Report();   
        
        tester.MeshNum++;
        tester.Run();
        EventHandler::Report();        
    }
};

#endif /*TESTPERFORMANCE_HPP_*/
