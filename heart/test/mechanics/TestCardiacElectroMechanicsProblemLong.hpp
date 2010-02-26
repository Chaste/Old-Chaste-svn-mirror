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


#ifndef TESTCARDIACELECTROMECHANICSPROBLEMLONG_HPP_
#define TESTCARDIACELECTROMECHANICSPROBLEMLONG_HPP_

#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include <petscvec.h>
#include "PetscSetupAndFinalize.hpp"
#include "CardiacElectroMechProbRegularGeom.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "NonlinearElasticityTools.hpp"

class TestCardiacElectroMechanicsProblemLong : public CxxTest::TestSuite
{
public:
    void Test2dHardcodedResult() throw(Exception)
    {
        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 2> cell_factory(-1000*1000);

        // run to 125 ms - about where the width is at its minimum (see figures
        // in "A numerical method for cardiac mechano–electric simulations" (Pras&JonW))
        CardiacElectroMechProbRegularGeom<2> problem(NHS,
                                                     1.0,  /* width */
                                                     5,    /* mech mesh size */
                                                     60,   /* elec elem each dir */
                                                     &cell_factory,
                                                     125,  /* end time */
                                                     0.01, /* electrics timestep (ms) */
                                                     100,  /* 100*0.01ms mech dt */
                                                     1.0,  /* contraction model ode dt */
                                                     "TestCardiacEmNhs2dLong");

        problem.SetNoElectricsOutput();
        problem.Solve();

        // test by checking the length of the tissue against hardcoded value
        std::vector<c_vector<double,2> >& r_deformed_position = problem.rGetDeformedPosition();
        TS_ASSERT_DELTA(r_deformed_position[5](0), 0.8282, 1e-3);

        MechanicsEventHandler::Headings();
        MechanicsEventHandler::Report();
    }

    void Test2dVariableFibres() throw(Exception)
    {
        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 2> cell_factory(-1000*1000);

        CardiacElectroMechProbRegularGeom<2> problem(NHS,
                                                     1.0,  /* width */
                                                     5,    /* mech mesh size */
                                                     96,   /* elec elem each dir */
                                                     &cell_factory,
                                                     125,  /* end time */
                                                     0.01, /* electrics timestep (ms) */
                                                     100,  /* 100*0.01ms mech dt */
                                                     1.0,  /* contraction model ode dt */
                                                     "TestCardiacEmVaryingFibres");
        
        // fibres going from (1,0) at X=0 to (1,1)-direction at X=1
        /* the fibres file was created with the code (inside a class that owns a mesh)
        for(unsigned elem_index=0; elem_index<this->mpQuadMesh->GetNumElements(); elem_index++)
        {
            double X = this->mpQuadMesh->GetElement(elem_index)->CalculateCentroid()[0];
            double theta = M_PI*X/4;
            std::cout << cos(theta) << " " << sin(theta) << " " << -sin(theta) << " " << cos(theta) << "\n" << std::flush;
        }
        assert(0);
        */
        problem.SetVariableFibreSheetDirectionsFile("heart/test/data/5by5mesh_curving_fibres.ortho");

        // problem.SetNoElectricsOutput();
        problem.Solve();

        // test by checking the length of the tissue against hardcoded value
        std::vector<c_vector<double,2> >& r_deformed_position = problem.rGetDeformedPosition();
        TS_ASSERT_DELTA(r_deformed_position[5](0), 0.9017, 2e-3); // visualised, looks good
        //IntelProduction differs by about 1.6e-3...
        
        MechanicsEventHandler::Headings();
        MechanicsEventHandler::Report();
    }



    void Test3d() throw(Exception)
    {
        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 3> cell_factory(-1000*1000);

        // set up two meshes of 1mm by 1mm by 1mm
        TetrahedralMesh<3,3> electrics_mesh;
        electrics_mesh.ConstructCuboid(10,10,10);
        electrics_mesh.Scale(0.01, 0.01, 0.01);

        QuadraticMesh<3> mechanics_mesh(0.1, 0.1, 0.1, 1, 1, 1);

        // fix the nodes on x=0
        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<3>::GetNodesByComponentValue(mechanics_mesh,0,0);

        CardiacElectroMechanicsProblem<3> problem(NHS,
                                                  &electrics_mesh,
                                                  &mechanics_mesh,
                                                  fixed_nodes,
                                                  &cell_factory,
                                                  50,   /* end time */
                                                  0.01, /* electrics timestep (ms) */
                                                  100,  /* 100*0.01ms mech dt */
                                                  1.0,  /* contraction model ode dt */
                                                  "TestCardiacElectroMech3d");

        problem.SetNoElectricsOutput();
        problem.Solve();

        // test by checking the length of the tissue against hardcoded value
        std::vector<c_vector<double,3> >& r_deformed_position = problem.rGetDeformedPosition();
        TS_ASSERT_DELTA(r_deformed_position[1](0), 0.0879, 1e-3);

        MechanicsEventHandler::Headings();
        MechanicsEventHandler::Report();
    }

};
#endif /*TESTCARDIACELECTROMECHANICSPROBLEMLONG_HPP_*/
