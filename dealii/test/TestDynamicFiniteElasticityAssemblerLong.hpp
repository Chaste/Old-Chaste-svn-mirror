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


#ifndef TESTDYNAMICFINITEELASTICITYASSEMBLERLONG_HPP_
#define TESTDYNAMICFINITEELASTICITYASSEMBLERLONG_HPP_


#include <cxxtest/TestSuite.h>
#include "DynamicFiniteElasticityAssembler.cpp"
#include "TriangulationVertexIterator.hpp"
#include "DofVertexIterator.hpp"
#include "MooneyRivlinMaterialLaw.hpp"
#include "PolynomialMaterialLaw3d.hpp"
#include "ExponentialMaterialLaw.hpp"
#include "FiniteElasticityTools.hpp"

// todos: proper test of answers

class TestDynamicFiniteElasticityAssemblerLong : public CxxTest::TestSuite
{
public :

    void Test2dProblemOnSquareLongerRun() throw(Exception)
    {
        Vector<double> body_force(2);
        body_force(1) = 2.0;

        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(2.0);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(3);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 0.0);


        DynamicFiniteElasticityAssembler<2> dynamic_finite_elasticity(&mesh,
                                                                      &mooney_rivlin_law,
                                                                      body_force,
                                                                      1.0,
                                                                      "dynamic_finite_elas/simple2dlong"
                                                                      );

        dynamic_finite_elasticity.SetTimes(0.0,1.0,0.01);

        dynamic_finite_elasticity.Solve();

        // get undeformed position
        std::vector<Vector<double> >& undeformed_position
            = dynamic_finite_elasticity.rGetUndeformedPosition();

        // get deformed position
        std::vector<Vector<double> >& deformed_position
            = dynamic_finite_elasticity.rGetDeformedPosition();

        for (unsigned vertex_index=0; vertex_index<deformed_position[0].size(); vertex_index++)
        {
            // todo: TEST THESE!!
            double X = undeformed_position[0](vertex_index);
            double Y = undeformed_position[1](vertex_index);
            double x = deformed_position[0](vertex_index);
            double y = deformed_position[1](vertex_index);
            std::cout << vertex_index << " " << X << " " << Y
                      << " " << x << " " << y << "\n";
        }
    }


    void Test3dProblemOnCube() throw(Exception)
    {
        Vector<double> body_force(3);
        body_force(1) = 0.02;

        MooneyRivlinMaterialLaw<3> mooney_rivlin_law(0.01, 0.02);

        Triangulation<3> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(3);
        FiniteElasticityTools<3>::SetFixedBoundary(mesh, 0, 0.0);


        DynamicFiniteElasticityAssembler<3> dynamic_finite_elasticity(&mesh,
                                                                      &mooney_rivlin_law,
                                                                      body_force,
                                                                      1.0,
                                                                      "dynamic_finite_elas/simple3d"
                                                                      );

        dynamic_finite_elasticity.SetTimes(0.0,0.2,0.01);

        dynamic_finite_elasticity.Solve();

        // get undeformed position
        std::vector<Vector<double> >& undeformed_position
            = dynamic_finite_elasticity.rGetUndeformedPosition();

        // get deformed position
        std::vector<Vector<double> >& deformed_position
            = dynamic_finite_elasticity.rGetDeformedPosition();

        for (unsigned vertex_index=0; vertex_index<deformed_position[0].size(); vertex_index++)
        {
            // todo: TEST THESE!!
            double X = undeformed_position[0](vertex_index);
            double Y = undeformed_position[1](vertex_index);
            double x = deformed_position[0](vertex_index);
            double y = deformed_position[1](vertex_index);
            std::cout << vertex_index << " " << X << " " << Y
                      << " " << x << " " << y << "\n";
        }
    }


    // Test that the final deformed mesh of a dynamic simulation which reaches a
    // resting state is the same as the result of a static simulation (ie the
    // final deformed mesh satisfies the static equations).
    void TestDynamicVsStatic2dOnSquare() throw(Exception)
    {
        Vector<double> body_force(2);
        body_force(1) = 2.0;
        body_force(1) = 1.0;

        double density = 1.0;
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(2.0);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(3);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 0.0);

        DynamicFiniteElasticityAssembler<2> dynamic_finite_elasticity(&mesh,
                                                                      &mooney_rivlin_law,
                                                                      body_force,
                                                                      density,
                                                                      "dynamic_finite_elas/test_dymamic_v_static"
                                                                     );

        dynamic_finite_elasticity.SetTimes(0.0,1.0,0.01);
        dynamic_finite_elasticity.Solve();

        FiniteElasticityAssembler<2> finite_elasticity(&mesh,
                                                       &mooney_rivlin_law,
                                                       body_force,
                                                       density,
                                                       "finite_elas/test_dymamic_v_static"
                                                      );
        finite_elasticity.StaticSolve();

        Vector<double>& dynamic_solution = dynamic_finite_elasticity.rGetCurrentSolution();
        Vector<double>& static_solution = finite_elasticity.rGetCurrentSolution();

        for (unsigned i=0; i<dynamic_solution.size(); i++)
        {
            double tol;

            // quick dirty check to see if solution(i) corresponds to
            // displacement or pressure:
            if (fabs(dynamic_solution(i))<1)
            {
                //probably displacement
                tol = 5e-3;
            }
            else
            {
                //probably pressure
                tol = 1e-1;
            }

            TS_ASSERT_DELTA(dynamic_solution(i), static_solution(i), tol);
        }
    }
};

#endif /*TESTDYNAMICFINITEELASTICITYASSEMBLERLONG_HPP_*/
