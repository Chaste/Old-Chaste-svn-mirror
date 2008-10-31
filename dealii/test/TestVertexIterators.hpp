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


#ifndef TESTVERTEXITERATORS_HPP_
#define TESTVERTEXITERATORS_HPP_

#include <cxxtest/TestSuite.h>

#include "FiniteElasticityAssembler.cpp"
#include "TriangulationVertexIterator.hpp"
#include "DofVertexIterator.hpp"
#include "MooneyRivlinMaterialLaw2.hpp"

#include "FiniteElasticityTools.hpp"

#include <grid/grid_generator.h>

class TestVertexIterators : public CxxTest::TestSuite
{
public :
    void TestTriangulationVertexIterator() throw(Exception)
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(1);

        TriangulationVertexIterator<2>  vertex_iter(&mesh);
        unsigned initial_index = vertex_iter.GetVertexGlobalIndex();

        unsigned counter=0;
        unsigned num_nodes = mesh.n_vertices();

        unsigned num_x_is_zero = 0;
        unsigned num_x_is_half = 0;
        unsigned num_x_is_one  = 0;

        unsigned num_y_is_zero = 0;
        unsigned num_y_is_half = 0;
        unsigned num_y_is_one  = 0;

        std::vector<bool> vertex_touched(mesh.n_vertices(),false);

        while (!vertex_iter.ReachedEnd())
        {
            Triangulation<2>::active_cell_iterator cell = vertex_iter.GetCell();

            //unsigned local_vertex_index = vertex_iter.GetLocalVertexIndexForCell();

            unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();

            TS_ASSERT(vertex_touched[vertex_index]==false);
            vertex_touched[vertex_index] = true;

            unsigned cell_vertex_index = vertex_iter.GetCell()->vertex_index(vertex_iter.GetLocalVertexIndexForCell());
            TS_ASSERT_EQUALS(cell_vertex_index, vertex_index);

            Point<2> vertex = vertex_iter.GetVertex();

            if (fabs(vertex(0)-0)<1e-6)
            {
                num_x_is_zero++;
            }
            else if (fabs(vertex(0)-0.5)<1e-6)
            {
                num_x_is_half++;
            }
            else if (fabs(vertex(0)-1)<1e-6)
            {
                num_x_is_one++;
            }

            if (fabs(vertex(1)-0)<1e-6)
            {
                num_y_is_zero++;
            }
            else if (fabs(vertex(1)-0.5)<1e-6)
            {
                num_y_is_half++;
            }
            else if (fabs(vertex(1)-1)<1e-6)
            {
                num_y_is_one++;
            }

            assert(counter++ < num_nodes);
            vertex_iter.Next();
        }

        TS_ASSERT_EQUALS(counter, num_nodes);

        TS_ASSERT_EQUALS(num_x_is_zero,3u);
        TS_ASSERT_EQUALS(num_x_is_half,3u);
        TS_ASSERT_EQUALS(num_x_is_one, 3u);
        TS_ASSERT_EQUALS(num_y_is_zero,3u);
        TS_ASSERT_EQUALS(num_y_is_half,3u);
        TS_ASSERT_EQUALS(num_y_is_one, 3u);

        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT(vertex_touched[i]==true);
        }

        vertex_iter.Reset();
        counter=0;
        TS_ASSERT_EQUALS(vertex_iter.GetVertexGlobalIndex(), initial_index);
        while (!vertex_iter.ReachedEnd())
        {
            counter++;
            vertex_iter.Next();
        }
        TS_ASSERT_EQUALS(counter, num_nodes);

        // do it again after refining the mesh, check no errors
        mesh.refine_global(1);
        vertex_iter.Reset();
        counter=0;
        TS_ASSERT_EQUALS(vertex_iter.GetVertexGlobalIndex(), initial_index);
        while (!vertex_iter.ReachedEnd())
        {
            counter++;
            vertex_iter.Next();
        }
        TS_ASSERT_EQUALS(counter, mesh.n_vertices());
    }




    void TestDofVertexIterator() throw(Exception)
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(1);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 0.0);

        // create a FiniteElasticity object, just in order to get its DoFHandler
        Vector<double> body_force(2);
        body_force(0) = 6.0;
        MooneyRivlinMaterialLaw2<2> mooney_rivlin_law(2.0,2.0);
        FiniteElasticityAssembler<2> finite_elasticity(&mesh, &mooney_rivlin_law,
                                                       body_force, 1.0, "finite_elas/simple");


        DoFHandler<2>& dof_handler = finite_elasticity.rGetDofHandler();

        DofVertexIterator<2>  vertex_iter(&mesh,&dof_handler);
        unsigned initial_index = vertex_iter.GetVertexGlobalIndex();

        unsigned counter=0;
        unsigned num_nodes = mesh.n_vertices();

        unsigned num_x_is_zero = 0;
        unsigned num_x_is_half = 0;
        unsigned num_x_is_one  = 0;

        unsigned num_y_is_zero = 0;
        unsigned num_y_is_half = 0;
        unsigned num_y_is_one  = 0;

        std::vector<bool> vertex_touched(mesh.n_vertices(),false);
        std::vector<bool> dof_touched(dof_handler.n_dofs(),false);

        while (!vertex_iter.ReachedEnd())
        {
            Triangulation<2>::active_cell_iterator cell = vertex_iter.GetCell();

            //unsigned local_vertex_index = vertex_iter.GetLocalVertexIndexForCell();

            unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
            TS_ASSERT(vertex_touched[vertex_index]==false);
            vertex_touched[vertex_index] = true;

            unsigned dof_index_for_x = vertex_iter.GetDof(0);
            unsigned dof_index_for_y = vertex_iter.GetDof(1);

            TS_ASSERT(dof_touched[dof_index_for_x]==false);
            TS_ASSERT(dof_touched[dof_index_for_y]==false);

            dof_touched[dof_index_for_x] = true;
            dof_touched[dof_index_for_y] = true;

            Point<2> vertex = vertex_iter.GetVertex();

            if (fabs(vertex(0)-0)<1e-6)
            {
                num_x_is_zero++;
            }
            else if (fabs(vertex(0)-0.5)<1e-6)
            {
                num_x_is_half++;
            }
            else if (fabs(vertex(0)-1)<1e-6)
            {
                num_x_is_one++;
            }

            if (fabs(vertex(1)-0)<1e-6)
            {
                num_y_is_zero++;
            }
            else if (fabs(vertex(1)-0.5)<1e-6)
            {
                num_y_is_half++;
            }
            else if (fabs(vertex(1)-1)<1e-6)
            {
                num_y_is_one++;
            }

            assert(counter++ < num_nodes);
            vertex_iter.Next();
        }

        TS_ASSERT_EQUALS(counter, num_nodes);

        TS_ASSERT_EQUALS(num_x_is_zero,3u);
        TS_ASSERT_EQUALS(num_x_is_half,3u);
        TS_ASSERT_EQUALS(num_x_is_one, 3u);
        TS_ASSERT_EQUALS(num_y_is_zero,3u);
        TS_ASSERT_EQUALS(num_y_is_half,3u);
        TS_ASSERT_EQUALS(num_y_is_one, 3u);

        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT(vertex_touched[i]==true);
        }

        vertex_iter.Reset();
        counter=0;
        TS_ASSERT_EQUALS(vertex_iter.GetVertexGlobalIndex(), initial_index);
        while (!vertex_iter.ReachedEnd())
        {
            counter++;
            vertex_iter.Next();
        }
        TS_ASSERT_EQUALS(counter, num_nodes);

        // do it again after refining the mesh, check no errors
        mesh.refine_global(1);
        vertex_iter.Reset();
        counter=0;
        TS_ASSERT_EQUALS(vertex_iter.GetVertexGlobalIndex(), initial_index);
        while (!vertex_iter.ReachedEnd())
        {
            counter++;
            vertex_iter.Next();
        }
        TS_ASSERT_EQUALS(counter, mesh.n_vertices());
    }
};

#endif /*TESTVERTEXITERATORS_HPP_*/
