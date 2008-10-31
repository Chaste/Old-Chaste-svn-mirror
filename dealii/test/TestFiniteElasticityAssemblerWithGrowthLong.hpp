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


#ifndef TESTFINITEELASTICITYASSEMBLERWITHGROWTHLONG_HPP_
#define TESTFINITEELASTICITYASSEMBLERWITHGROWTHLONG_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include "FiniteElasticityAssemblerWithGrowth.cpp"
#include "TriangulationVertexIterator.hpp"
#include "MooneyRivlinMaterialLaw2.hpp"
#include "FiniteElasticityTools.hpp"
#include "ConstantTumourSourceModel.hpp"



class TestFiniteElasticityAssemblerWithGrowthLong : public CxxTest::TestSuite
{
public :

    void TestWithSimpleProblem3d() throw(Exception)
    {
        Vector<double> body_force(3); // zero
        double density = 1.233;

        MooneyRivlinMaterialLaw2<3> mooney_rivlin_law(0.02, 0.01);

        Triangulation<3> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(2);

        double initial_elem_volume = 1.0/mesh.n_active_cells();

        Point<3> zero;
        FiniteElasticityTools<3>::FixFacesContainingPoint(mesh, zero);

        // set all elements as growing
        FiniteElasticityTools<3>::SetCircularRegionAsGrowingRegion(mesh, zero, 100);
        double source_value = 2;
        ConstantTumourSourceModel<3> source_model(source_value);

        FiniteElasticityAssemblerWithGrowth<3> finiteelas_with_growth(&mesh,
                                                                      &mooney_rivlin_law,
                                                                      body_force,
                                                                      density,
                                                                      "finite_elas_growth/simple3d",
                                                                      &source_model);


        // loop over all the elements, and if it is in the growing region, check
        // each node has an ode system associated with it...
        TriangulationVertexIterator<3> vertex_iter(&mesh);
        while(!vertex_iter.ReachedEnd())
        {
            unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
            TS_ASSERT_EQUALS(finiteelas_with_growth.IsGrowingNode(vertex_index), true);
            vertex_iter.Next();
        }

        // run
        double end_time = 0.2;
        finiteelas_with_growth.SetTimes(0.0, end_time, 0.1);
        finiteelas_with_growth.Run();

        std::vector<Vector<double> >& deformed_position
            = finiteelas_with_growth.rGetDeformedPosition();

        // test
        Triangulation<3>::active_cell_iterator element_iter = mesh.begin_active();
        while (element_iter!=mesh.end())
        {
            double x_undef_node_0 = element_iter->vertex(0)[0];
            double y_undef_node_0 = element_iter->vertex(0)[1];
            double z_undef_node_0 = element_iter->vertex(0)[2];

            // look at elements in the quadrant containing (1,1,1), ie those away
            // from the fixed corner (0,0,0). (if there was no fixed nodes, then
            // the solution would just be iostropic enlargement (simple stretching)
            if( (x_undef_node_0>0.5) && (y_undef_node_0>0.5) && (z_undef_node_0>0.5) )
            {
                std::vector<std::vector<double> > elem_def_posn(3);
                elem_def_posn[0].resize(8);
                elem_def_posn[1].resize(8);
                elem_def_posn[2].resize(8);

                for(unsigned i=0; i<3; i++)
                {
                    for(unsigned j=0; j<8; j++)
                    {
                        elem_def_posn[i][j] = deformed_position[i](element_iter->vertex_index(j));
                    }
                }

                /* Cube node ordering
                 * node 0: (0 0 0)
                 * node 1: (1 0 0)
                 * node 2: (1 0 1)
                 * node 3: (0 0 1)
                 * node 4: (0 1 0)
                 * node 5: (1 1 0)
                 * node 6: (1 1 1)
                 * node 7: (0 1 1)
                 */

                // these elements are away from the fixed corner, so should be very
                // like enlarged rectangles. Verify this
                // x0 = x3 = x4 = x7
                TS_ASSERT_DELTA(elem_def_posn[0][0], elem_def_posn[0][3], 1e-3);
                TS_ASSERT_DELTA(elem_def_posn[0][0], elem_def_posn[0][4], 1e-3);
                TS_ASSERT_DELTA(elem_def_posn[0][0], elem_def_posn[0][7], 1e-3);

                // x1 = x2 = x5 = x6
                TS_ASSERT_DELTA(elem_def_posn[0][1], elem_def_posn[0][2], 1e-3);
                TS_ASSERT_DELTA(elem_def_posn[0][1], elem_def_posn[0][5], 1e-3);
                TS_ASSERT_DELTA(elem_def_posn[0][1], elem_def_posn[0][6], 1e-3);

                // y0 = y1 = y2 = y3
                TS_ASSERT_DELTA(elem_def_posn[1][0], elem_def_posn[1][1], 1e-3);
                TS_ASSERT_DELTA(elem_def_posn[1][0], elem_def_posn[1][2], 1e-3);
                TS_ASSERT_DELTA(elem_def_posn[1][0], elem_def_posn[1][3], 1e-3);

                // y4 = y5 = y6 = y7
                TS_ASSERT_DELTA(elem_def_posn[1][4], elem_def_posn[1][5], 1e-3);
                TS_ASSERT_DELTA(elem_def_posn[1][4], elem_def_posn[1][6], 1e-3);
                TS_ASSERT_DELTA(elem_def_posn[1][4], elem_def_posn[1][7], 1e-3);

                // z0 = z1 = z4 = z5
                TS_ASSERT_DELTA(elem_def_posn[2][0], elem_def_posn[2][1], 1e-3);
                TS_ASSERT_DELTA(elem_def_posn[2][0], elem_def_posn[2][4], 1e-3);
                TS_ASSERT_DELTA(elem_def_posn[2][0], elem_def_posn[2][5], 1e-3);

                // z2 = z3 = z6 = z7
                TS_ASSERT_DELTA(elem_def_posn[2][2], elem_def_posn[2][3], 1e-3);
                TS_ASSERT_DELTA(elem_def_posn[2][2], elem_def_posn[2][6], 1e-3);
                TS_ASSERT_DELTA(elem_def_posn[2][2], elem_def_posn[2][7], 1e-3);


                // check volume of enlarged square is as expected
                // dgdt = 1/3 g rho s, so g = exp(1/3 rho s T)
                // F = gI so detF = g^3 = exp(rho s T)
                double expected_volume = exp(density*end_time*source_value)*initial_elem_volume;
                double volume =   (elem_def_posn[0][1]-elem_def_posn[0][0])  // x1-x0
                                * (elem_def_posn[1][4]-elem_def_posn[1][0])  // y4-y0
                                * (elem_def_posn[2][3]-elem_def_posn[2][0]); // z3-z0
                TS_ASSERT_DELTA( volume, expected_volume, 1e-3);
            }
            element_iter++;
        }
    }


    void TestWithSimpleProblemForRefinement() throw(Exception)
    {
        Vector<double> body_force(2); // zero
        double density = 1;

        MooneyRivlinMaterialLaw2<2> mooney_rivlin_law(0.02);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(3);

        // for a second simulation without refinement
        Triangulation<2> mesh2;
        GridGenerator::hyper_cube(mesh2, 0.0, 1.0);
        mesh2.refine_global(3);

        double initial_elem_volume = 1.0/mesh.n_active_cells();

        Point<2> zero;
        FiniteElasticityTools<2>::FixFacesContainingPoint(mesh, zero);
        FiniteElasticityTools<2>::FixFacesContainingPoint(mesh2, zero);

        // set all elements as growing
        FiniteElasticityTools<2>::SetCircularRegionAsGrowingRegion(mesh, zero, 100);
        FiniteElasticityTools<2>::SetCircularRegionAsGrowingRegion(mesh2, zero, 100);

        double source_value = 2;
        ConstantTumourSourceModel<2> source_model(source_value);
        ConstantTumourSourceModel<2> source_model2(source_value);

        FiniteElasticityAssemblerWithGrowth<2> finiteelas_with_growth(&mesh,
                                                                      &mooney_rivlin_law,
                                                                      body_force,
                                                                      density,
                                                                      "finite_elas_growth/simple2drefinement",
                                                                      &source_model);
        // a second one that does not remesh
        FiniteElasticityAssemblerWithGrowth<2> finiteelas_with_growth2(&mesh2,
                                                                       &mooney_rivlin_law,
                                                                       body_force,
                                                                       density,
                                                                       "finite_elas_growth/simple2drefinement2",
                                                                       &source_model2);
        finiteelas_with_growth2.SetNoRefinement();

        // run
        double end_time = 0.8;
        finiteelas_with_growth.SetTimes(0.0, end_time, 0.1);
        finiteelas_with_growth.Run();

        finiteelas_with_growth2.SetTimes(0.0, end_time, 0.1);
        finiteelas_with_growth2.Run();


        std::vector<Vector<double> >& deformed_position
            = finiteelas_with_growth.rGetDeformedPosition();

        // test
        Triangulation<2>::active_cell_iterator element_iter = mesh.begin_active();
        while (element_iter!=mesh.end())
        {
            double x_undef_node_0 = element_iter->vertex(0)[0];
            double y_undef_node_0 = element_iter->vertex(0)[1];

            // look at elements in the top right-hand quadrant, ie those away
            // from the fixed corner (0,0). (if there was no fixed nodes, then
            // the solution would just be iostropic enlargement (simple stretching)
            if( (x_undef_node_0>0.5) && (y_undef_node_0>0.5))
            {
                double x0 = deformed_position[0](element_iter->vertex_index(0));
                double y0 = deformed_position[1](element_iter->vertex_index(0));

                double x1 = deformed_position[0](element_iter->vertex_index(1));
                double y1 = deformed_position[1](element_iter->vertex_index(1));

                double x2 = deformed_position[0](element_iter->vertex_index(2));
                double y2 = deformed_position[1](element_iter->vertex_index(2));

                double x3 = deformed_position[0](element_iter->vertex_index(3));
                double y3 = deformed_position[1](element_iter->vertex_index(3));

                // these elements are away from the fixed corner, so should be very
                // like enlarged rectangles. Verify this, ie check x0=x3, x1=x2, etc
                TS_ASSERT_DELTA(x0, x3, 1e-3);
                TS_ASSERT_DELTA(x1, x2, 1e-3);
                TS_ASSERT_DELTA(y0, y1, 1e-3);
                TS_ASSERT_DELTA(y2, y3, 1e-3);

                // check volume of enlarged square is as expected
                // dgdt = 1/2 g rho s, so g = exp(1/2 rho s T)
                // F = gI so detF = g^2 = exp(rho s T)
                double expected_volume = exp(density*end_time*source_value)*initial_elem_volume;

                // make sure we ran long enough for there to be refinement
                TS_ASSERT_LESS_THAN(4.0, expected_volume/initial_elem_volume)

                double volume = (x1-x0)*(y3-y0);

                // divide by 4 as there will have been refinement
                // Note - we run for quite long in this simulation (in order to get
                // refinement), so error in ode solving/fixed corner effects means the
                // volumes near the other corner are not as close as might be liked
                // (0.0101 != 0.0126). Note there is an extra test below
                TS_ASSERT_DELTA( volume, expected_volume/4, 5e-3);
            }
            element_iter++;
        }

        // compare the solution at (1,1) which the solution of the simulation without
        // refinement
        std::vector<Vector<double> >& deformed_position2
            = finiteelas_with_growth2.rGetDeformedPosition();
        TS_ASSERT_DELTA(deformed_position[0](2), deformed_position2[0](2), 1e-2);
        TS_ASSERT_DELTA(deformed_position[1](2), deformed_position2[1](2), 1e-2);
    }

    void TestWithSimpleProblemForCoarsening() throw(Exception)
    {
        Vector<double> body_force(2); // zero
        double density = 1;

        MooneyRivlinMaterialLaw2<2> mooney_rivlin_law(0.02);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(3);

        // for a second simulation without coarsening
        Triangulation<2> mesh2;
        GridGenerator::hyper_cube(mesh2, 0.0, 1.0);
        mesh2.refine_global(3);

        double initial_elem_volume = 1.0/mesh.n_active_cells();

        Point<2> zero;
        FiniteElasticityTools<2>::FixFacesContainingPoint(mesh, zero);
        FiniteElasticityTools<2>::FixFacesContainingPoint(mesh2, zero);

        // set all elements as growing
        FiniteElasticityTools<2>::SetCircularRegionAsGrowingRegion(mesh, zero, 100);
        FiniteElasticityTools<2>::SetCircularRegionAsGrowingRegion(mesh2, zero, 100);

        double source_value = -2;
        ConstantTumourSourceModel<2> source_model(source_value);
        ConstantTumourSourceModel<2> source_model2(source_value);

        FiniteElasticityAssemblerWithGrowth<2> finiteelas_with_growth(&mesh,
                                                                      &mooney_rivlin_law,
                                                                      body_force,
                                                                      density,
                                                                      "finite_elas_growth/simple2dcoarsen",
                                                                      &source_model);

        FiniteElasticityAssemblerWithGrowth<2> finiteelas_with_growth2(&mesh2,
                                                                       &mooney_rivlin_law,
                                                                       body_force,
                                                                       density,
                                                                       "finite_elas_growth/simple2dcoarsen2",
                                                                       &source_model2);


        // run
        double end_time = 0.8;
        finiteelas_with_growth.SetTimes(0.0, end_time, 0.1);
        finiteelas_with_growth.Run();

        finiteelas_with_growth2.SetTimes(0.0, end_time, 0.1);
        finiteelas_with_growth2.Run();

        std::vector<Vector<double> >& deformed_position
            = finiteelas_with_growth.rGetDeformedPosition();

        // test
        Triangulation<2>::active_cell_iterator element_iter = mesh.begin_active();
        while (element_iter!=mesh.end())
        {
            double x_undef_node_0 = element_iter->vertex(0)[0];
            double y_undef_node_0 = element_iter->vertex(0)[1];

            // look at elements in the top right-hand quadrant, ie those away
            // from the fixed corner (0,0). (if there was no fixed nodes, then
            // the solution would just be iostropic enlargement (simple stretching)
            if( (x_undef_node_0>0.5) && (y_undef_node_0>0.5))
            {
                double x0 = deformed_position[0](element_iter->vertex_index(0));
                double y0 = deformed_position[1](element_iter->vertex_index(0));

                double x1 = deformed_position[0](element_iter->vertex_index(1));
                double y1 = deformed_position[1](element_iter->vertex_index(1));

                double x2 = deformed_position[0](element_iter->vertex_index(2));
                double y2 = deformed_position[1](element_iter->vertex_index(2));

                double x3 = deformed_position[0](element_iter->vertex_index(3));
                double y3 = deformed_position[1](element_iter->vertex_index(3));

                // these elements are away from the fixed corner, so should be very
                // like enlarged rectangles. Verify this, ie check x0=x3, x1=x2, etc
                TS_ASSERT_DELTA(x0, x3, 1e-3);
                TS_ASSERT_DELTA(x1, x2, 1e-3);
                TS_ASSERT_DELTA(y0, y1, 1e-3);
                TS_ASSERT_DELTA(y2, y3, 1e-3);

                // check volume of enlarged square is as expected
                // dgdt = 1/2 g rho s, so g = exp(1/2 rho s T)
                // F = gI so detF = g^2 = exp(rho s T)
                double expected_volume = exp(density*end_time*source_value)*initial_elem_volume;

                // make sure we ran long enough for there to be coarsening
                TS_ASSERT_LESS_THAN(expected_volume/initial_elem_volume, 1.0/4.0);

                double volume = (x1-x0)*(y3-y0);

                // multiply by 4 as there will have been coarsening
                // Note - we run for quite long in this simulation (in order to get
                // coarsening), so error in ode solving/fixed corner effects means the
                // volumes near the other corner are not as close as might be liked
                // (0.0101 != 0.0126). Note there is an extra test below
                TS_ASSERT_DELTA( volume, 4*expected_volume, 3e-3);
            }
            element_iter++;
        }


        // compare the solution at (1,1) which the solution of the simulation without
        // coarsening
        std::vector<Vector<double> >& deformed_position2
            = finiteelas_with_growth2.rGetDeformedPosition();
        TS_ASSERT_DELTA(deformed_position[0](2), deformed_position2[0](2), 1e-2);
        TS_ASSERT_DELTA(deformed_position[1](2), deformed_position2[1](2), 1e-2);
    }
};


#endif /*TESTFINITEELASTICITYASSEMBLERWITHGROWTHLONG_HPP_*/


