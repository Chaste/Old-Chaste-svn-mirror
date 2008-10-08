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


#ifndef TESTCARDIACELECTROMECHANICSPROBLEM2_HPP_
#define TESTCARDIACELECTROMECHANICSPROBLEM2_HPP_

#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include <petscvec.h>
#include "PetscSetupAndFinalize.hpp"
#include "CardiacElectroMechanicsProblem2.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"

class TestCardiacElectroMechanicsProblem2 : public CxxTest::TestSuite
{
public:

    void TestDeterminingWatchedNodes() throw(Exception)
    {
        EventHandler::Disable();

        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 2> cell_factory(-1000*1000);

        CardiacElectroMechanicsProblem2<2> problem(&cell_factory,
                                                   1, /* end time */
                                                   1, /*mech mesh size*/
                                                   100, /* 100*0.01ms mech dt */
                                                   0.01,
                                                   "nothingtolookathere");

        c_vector<double,2> pos;
        pos(0) = 1.0;
        pos(1) = 0.0;
        problem.SetWatchedPosition(pos);
        problem.Initialise();

        // have checked these hardcoded values correspond to the nodes
        // at (1,0);
        TS_ASSERT_EQUALS(problem.mWatchedElectricsNodeIndex, 9408u);
        TS_ASSERT_EQUALS(problem.mWatchedMechanicsNodeIndex, 1u);

        //// would like to do the following....
        //CardiacElectroMechanicsProblem<2> problem2(&cell_factory,
        //                                           1, 10, 100, 0.01,
        //                                           "nothingtolookathere");
        //pos(1) = 1.1;
        //problem2.SetWatchedPosition(pos);
        //TS_ASSERT_THROWS_ANYTHING(problem2.Initialise());
        //// ... but the exception causes a segmentation fault and had to be replaced
        //// with an assert(0);
    }


    void Test2dOneMechanicsElement() throw(Exception)
    {
        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 2> cell_factory(-1000*1000);

        CardiacElectroMechanicsProblem2<2> problem(&cell_factory,
                                                   10, /* end time */
                                                   1, /*mech mesh size*/
                                                   100, /* 100*0.01ms mech dt */
                                                   0.01, /*NHS ode timestep */
                                                   "TestCardiacElectroMechOneElement",
                                                   0.05, /* Width of tissue */
                                                   5);   /* Num electrics elements in each dir */
        problem.SetNoElectricsOutput();

        c_vector<double,2> pos;
        pos(0) = 0.05;
        pos(1) = 0.0;
        
        problem.SetWatchedPosition(pos);
        problem.Solve();
 
        // test by checking the length of the tissue against hardcoded value
        std::vector<c_vector<double,2> >& r_deformed_position = problem.rGetDeformedPosition();
        TS_ASSERT_DELTA(r_deformed_position[1](0), 0.0497, 1e-4);

        OutputFileHandler handler("TestCardiacElectroMechOneElement",false);
        std::string watched = handler.GetOutputDirectoryFullPath() + "watched.txt";
        std::string command = "diff " + handler.GetOutputDirectoryFullPath() + "watched.txt heart/test/data/good_watched.txt";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);
        
        // check no electrics output was written
        command = "ls " + handler.GetOutputDirectoryFullPath() + "/electrics";
        TS_ASSERT_DIFFERS(system(command.c_str()), 0);
    }

//    void TestCinverseDataStructure() throw(Exception)
//    {
//        PlaneStimulusCellFactory<2> cell_factory(0.01, -1000*1000);
//        CardiacElectroMechanicsProblem<2> implicit_problem(&cell_factory,
//                                                           0.05, /* end time */
//                                                           5, /*mech mesh size*/
//                                                           1, /* 0.01ms mech dt */
//                                                           0.01,
//                                                           "TestCardiacElectroMechImplicitCinverse");
//        implicit_problem.SetNoElectricsOutput();
//        implicit_problem.Solve();
//
//        NodewiseData<2>* p_nodewise_data = NodewiseData<2>::Instance();
//        unsigned num_nodes = 97*97; //hardcoded
//        TS_ASSERT_EQUALS(p_nodewise_data->rGetData().size(), num_nodes);
//
//        for(unsigned i=0; i<p_nodewise_data->rGetData().size(); i++)
//        {
//            TS_ASSERT_EQUALS(p_nodewise_data->rGetData()[i].size(), 3u);
//            TS_ASSERT_DELTA(p_nodewise_data->rGetData()[i][0], 1.0, 1e-6);
//            TS_ASSERT_DELTA(p_nodewise_data->rGetData()[i][1], 0.0, 1e-6);
//            TS_ASSERT_DELTA(p_nodewise_data->rGetData()[i][2], 1.0, 1e-6);
//        }
//    }
};
#endif /*TESTCARDIACELECTROMECHANICSPROBLEM2_HPP_*/
