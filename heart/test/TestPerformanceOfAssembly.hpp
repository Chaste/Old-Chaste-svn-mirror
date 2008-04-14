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

#ifndef TESTPERFORMANCEOFASSEMBLY_HPP_
#define TESTPERFORMANCEOFASSEMBLY_HPP_

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

class TestPerformanceOfAssembly : public CxxTest::TestSuite
{   
public:

    
    void TestPerf() throw(Exception)
    {
        // write headings
        PerformanceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<2>, 2>::DisplayHeadings();
        EventHandler::Headings();

        // base line
        PerformanceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<2>, 2> tester;
        tester.SimTime=0.0025;
        
        for (unsigned mesh_num=0; mesh_num<3; mesh_num++)
        {
            tester.MeshNum=mesh_num;
            tester.Run();
            EventHandler::Report();      
        }
        
    }
};

#endif /*TESTPERFORMANCE_HPP_*/
