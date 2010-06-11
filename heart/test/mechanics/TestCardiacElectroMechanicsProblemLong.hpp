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
#include "NobleVargheseKohlNoble1998WithSac.hpp"

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
        TS_ASSERT_DELTA(r_deformed_position[5](0), 0.8257, 1e-3);

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
        problem.SetVariableFibreSheetDirectionsFile("heart/test/data/5by5mesh_curving_fibres.ortho", false);

        // problem.SetNoElectricsOutput();
        problem.Solve();

        // test by checking the length of the tissue against hardcoded value
        std::vector<c_vector<double,2> >& r_deformed_position = problem.rGetDeformedPosition();
        // visualised, looks good - contracts in X-direction near the fixed surface,
        // but on the other side the fibres are in the (1,1) direction, so contraction
        // pulls the tissue downward a bit
        TS_ASSERT_DELTA(r_deformed_position[5](0), 0.8677, 2e-3); 
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
    
    void TestTwistingCube() throw(Exception)
    {
        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 3> cell_factory(-1000*1000);

        // set up two meshes of 1mm by 1mm by 1mm
        TetrahedralMesh<3,3> electrics_mesh;
        electrics_mesh.ConstructCuboid(10,10,10);
        electrics_mesh.Scale(0.01, 0.01, 0.01);

        QuadraticMesh<3> mechanics_mesh(0.1, 0.1, 0.1, 5, 5, 5);

        // fix the nodes on Z=0
        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<3>::GetNodesByComponentValue(mechanics_mesh,2,0.0);

        CardiacElectroMechanicsProblem<3> problem(KERCHOFFS2003,
                                                  &electrics_mesh,
                                                  &mechanics_mesh,
                                                  fixed_nodes,
                                                  &cell_factory,
                                                  50,   /* end time */
                                                  0.01, /* electrics timestep (ms) */
                                                  100,  /* 100*0.01ms mech dt */
                                                  1.0,  /* contraction model ode dt */
                                                  "TestCardiacElectroMech3dTwistingCube");


/////// Use the following to set up the fibres file        
////        GaussianQuadratureRule<3> quad_rule(3);
////        QuadraturePointsGroup<3> quad_points(mechanics_mesh, quad_rule);
////        std::cout << quad_points.Size() << "\n"; 
////        for(unsigned i=0; i<quad_points.Size(); i++)
////        {
////            ////std::cout << quad_points.Get(i)(0) << " " << quad_points.Get(i)(1) << " " << quad_points.Get(i)(2) << " ";
////            double x = quad_points.Get(i)(0);
////            double theta = M_PI/3 - 10*x*2*M_PI/3; // 60 degrees when x=0, -60 when x=0.1; 
////            std::cout <<  "0 " << cos(theta)  << " " << sin(theta) 
////                      << " 0 " << -sin(theta) << " " << cos(theta)
////                      << " 1 0 0\n";  
////        }

        problem.SetVariableFibreSheetDirectionsFile("heart/test/data/5by5by5_fibres_by_quadpt.orthoquad", true);

        problem.Solve();

        // verified that it twists by visualising, some hardcoded values here..

        std::vector<c_vector<double,3> >& r_deformed_position = problem.rGetDeformedPosition();
        TS_ASSERT_DELTA(r_deformed_position[6*6*5](0),  0.0195, 1e-3);
        TS_ASSERT_DELTA(r_deformed_position[6*6*5](1), -0.0190, 1e-3);
        TS_ASSERT_DELTA(r_deformed_position[6*6*5](2),  0.1021, 1e-3);

        TS_ASSERT_DELTA(r_deformed_position[6*6*6-1](0), 0.0812, 1e-3);
        TS_ASSERT_DELTA(r_deformed_position[6*6*6-1](1), 0.1150, 1e-3);
        TS_ASSERT_DELTA(r_deformed_position[6*6*6-1](2), 0.1030, 1e-3);

        MechanicsEventHandler::Headings();
        MechanicsEventHandler::Report();

    }



//    void Test3dWithNoble98SacAndImpact() throw(Exception)
//    {
//        // zero stimuli
//        PlaneStimulusCellFactory<CML_noble_varghese_kohl_noble_1998_basic_with_sac, 3> cell_factory(0);
//
//        // set up two meshes of 1cm by 1cm by 1cm
//        TetrahedralMesh<3,3> electrics_mesh;
//        unsigned num_elem = 10;
//        electrics_mesh.ConstructCuboid(num_elem,num_elem,num_elem);
//        electrics_mesh.Scale(1.0/num_elem, 1.0/num_elem, 1.0/num_elem);
//
//        QuadraticMesh<3> mechanics_mesh(1.0, 1.0, 1.0, 5, 5, 5);
//
//        // fix the nodes on x=0
//        std::vector<unsigned> fixed_nodes
//          = NonlinearElasticityTools<3>::GetNodesByComponentValue(mechanics_mesh,0,0);
//
//        std::vector<BoundaryElement<2,3>*> impact_region;        
//        for (TetrahedralMesh<3,3>::BoundaryElementIterator iter
//              = mechanics_mesh.GetBoundaryElementIteratorBegin();
//             iter != mechanics_mesh.GetBoundaryElementIteratorEnd();
//             ++iter)
//        {
//            c_vector<double,3> centroid =(*iter)->CalculateCentroid(); 
//            if (    (fabs(centroid[1])<1e-4) 
//                 && (centroid[0] < 0.05)
//                 && (centroid[2] < 0.05) )
//            {
//                BoundaryElement<2,3>* p_element = *iter;
//                impact_region.push_back(p_element);
//            }
//        }
//        assert(impact_region.size() > 0);
//
//        CardiacElectroMechanicsProblem<3> problem(KERCHOFFS2003,
//                                                  &electrics_mesh,
//                                                  &mechanics_mesh,
//                                                  fixed_nodes,
//                                                  &cell_factory,
//                                                  10,   /* end time */
//                                                  0.01, /* electrics timestep (ms) */
//                                                  100,  /* 100*0.01ms mech dt */
//                                                  1.0,  /* contraction model ode dt */
//                                                  "TestCardiacElectroMech3dImpact");
//
//        problem.SetImpactRegion(impact_region);
//
//        problem.Solve();
//
//        // test by checking the length of the tissue against hardcoded value
//        std::vector<c_vector<double,3> >& r_deformed_position = problem.rGetDeformedPosition();
//        TS_ASSERT_DELTA(r_deformed_position[1](0), 0.0879, 1e-3);
//
//        MechanicsEventHandler::Headings();
//        MechanicsEventHandler::Report();
//    }
};
#endif /*TESTCARDIACELECTROMECHANICSPROBLEMLONG_HPP_*/
