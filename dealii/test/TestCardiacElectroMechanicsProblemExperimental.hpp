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


#ifndef TESTCARDIACELECTROMECHANICSPROBLEMEXPERIMENTAL_HPP_
#define TESTCARDIACELECTROMECHANICSPROBLEMEXPERIMENTAL_HPP_

#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include <petscvec.h>
#include "PetscSetupAndFinalize.hpp"
#include "CardiacElectroMechanicsProblem1d.hpp"
#include "CardiacElectroMechanicsProblem.hpp"
#include "RegularStimulus.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"


class MyCellFactory : public AbstractCardiacCellFactory<2>
{
private:
    RegularStimulus* mpStimulus;

public:
    MyCellFactory() : AbstractCardiacCellFactory<2>()
    {
        mpStimulus = new RegularStimulus(-600, 0.5, 100000, 0.0);
        LOG(2, "Using RegularStimulus(-600, 0.5, 100000, 0.0)\n");
    }

    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        if (this->mpMesh->GetNode(node)->GetPoint()[0] == 0.0)
        {
            return new LuoRudyIModel1991OdeSystem(this->mpSolver,
                                                  mpStimulus,
                                                  this->mpZeroStimulus);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(this->mpSolver,
                                                  this->mpZeroStimulus,
                                                  this->mpZeroStimulus);
        }
    }

    ~MyCellFactory()
    {
        delete mpStimulus;
    }
};


class TestCardiacElectroMechanicsProblemExperimental : public CxxTest::TestSuite
{
public:

    void Test2dBasicRun() throw(Exception)
    {
        //PlaneStimulusCellFactory<2> cell_factory(0.01, -1000*1000);
        MyCellFactory cell_factory;

        CardiacElectroMechanicsProblem<2> implicit_problem(&cell_factory,
                                                           2300,  // end time of 200ms
                                                           10,   // 10 mech elements in 1cm by 1cm square
                                                           100,  // 100 mech times per elec timestep, ie mech_dt = 1ms
                                                           1,    // nhs_dt = 1ms
                                                           "CardiacElectroMechBasicLong");
        implicit_problem.SetNoElectricsOutput();
        implicit_problem.UseDirectLinearSolver();

        c_vector<double,2> pos;
        pos(0) = 1.0;
        pos(1) = 0.0;
        implicit_problem.SetWatchedPosition(pos);

        implicit_problem.Solve();
    }


//    void dontTestScaleCalcium() throw(Exception)
//    {
//        EventHandler::Disable();
//
//        double calcium_scale_factors[8] = {0.9, 0.95, 0.99, 1.0, 1.01, 1.05, 1.1, 2};
//
//        for(unsigned i=0; i<8; i++)
//        {
//            std::stringstream name;
//            name << "CardiacElectroMechScaleCalcium/" << calcium_scale_factors[i];
//
//            PlaneStimulusCellFactory<2> cell_factory(0.01, -1000*1000);
//
//            CardiacElectroMechanicsProblem<2> implicit_problem(&cell_factory,
//                                                               200,  // end time of 200ms
//                                                               10,   // 10 mech elements in 1cm by 1cm square
//                                                               100,  // 100 mech times per elec timestep, ie mech_dt = 1ms
//                                                               1,    // nhs_dt = 1ms
//                                                               name.str(),
//                                                               calcium_scale_factors[i]);
//            implicit_problem.SetNoElectricsOutput();
//            implicit_problem.Solve();
//        }
//    }



//    void Test2dCompareExplicitVsImplicit() throw(Exception)
//    {
//        PlaneStimulusCellFactory<2> cell_factory(0.01, -1000*1000);
//
//        unsigned dt=128;
//        for(unsigned i=0;i<4;i++)
//        {
//            std::stringstream name;
//            name << "CardiacElectroMech_Time_MORE_" << dt;
//
//            CardiacElectroMechanicsProblem<2> implicit_problem(&cell_factory, 200, 10, dt, 0.01, name.str());
//            implicit_problem.SetNoElectricsOutput();
//            implicit_problem.Solve();
//
//            dt *= 2;
//        }
//
//
//        unsigned num_nodes_in_each_dir=8;
//        for(unsigned i=0;i<1;i++)
//        {
//            std::stringstream name;
//            name << "CardiacElectroMech_Space_NEW_" << num_nodes_in_each_dir;
//
//            CardiacElectroMechanicsProblem<2> implicit_problem(&cell_factory, 200, num_nodes_in_each_dir, 100, 1.0, name.str());
//            implicit_problem.SetNoElectricsOutput();
//            implicit_problem.Solve();
//
//            num_nodes_in_each_dir *= 2;
//        }

//        double nhs_ode_time_step = 0.01;
//
//        for(unsigned i=0; i<8; i++)
//        {
//            std::stringstream name;
//            name << "CardiacElectroMech_OdeTimeStepNew_" << i;
//
//            std::cout << nhs_ode_time_step << "\n";
//
//            CardiacElectroMechanicsProblem<2> implicit_problem(&cell_factory, 200, 10, 128, nhs_ode_time_step, name.str());
//            implicit_problem.SetNoElectricsOutput();
//            implicit_problem.Solve();
//
//            nhs_ode_time_step *= 2;
//        }
//    }


};
#endif /*TESTCARDIACELECTROMECHANICSPROBLEMEXPERIMENTAL_HPP_*/
