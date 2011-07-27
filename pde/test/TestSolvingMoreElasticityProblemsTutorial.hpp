/*

Copyright (C) University of Oxford, 2005-2011

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

/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 */
#ifndef TESTSOLVINGMOREELASTICITYPROBLEMSTUTORIALS_HPP_
#define TESTSOLVINGMOREELASTICITYPROBLEMSTUTORIALS_HPP_

/* In this second solid mechanics tutorial, we illustrate some other possibilities. Specifically, we will
 * show the (very minor) changes required to solve a compressible nonlinear elasticity problem, we will
 * describe and show how to specify 'pressure on deformed body' boundary conditions, we illustrate how
 * a quadratic mesh can be generated using a linear mesh input files, and we also illustrate how
 * `Solve()` can be called repeatedly, with loading changing between the solves.
 *
 * EMPTYLINE
 *
 * These includes are the same as before
 */
#include <cxxtest/TestSuite.h>
#include "TrianglesMeshReader.hpp"
#include "PetscSetupAndFinalize.hpp"
/* These two are specific to compressible problems */
#include "CompressibleNonlinearElasticitySolver.hpp"
#include "CompressibleMooneyRivlinMaterialLaw.hpp"

class TestSolvingMoreElasticityProblemsTutorials : public CxxTest::TestSuite
{
public:
    void TestSolvingCompressibleProblem() throw (Exception)
    {
        /* All mechanics problems must take in quadratic meshes, but the mesh files for
         * (standard) linear meshes in Triangles/Tetgen can be automatically converted
         * to quadratic meshes, by simply doing the following. (The mesh loaded here is a disk
         * centred at the origin with radius 1).
         */
        QuadraticMesh<2> mesh;
        TrianglesMeshReader<2,2> reader("mesh/test/data/disk_522_elements");
        mesh.ConstructFromLinearMeshReader(reader);

        /* Compressible problems require a compressible material law, ie one that
         * inherits from `AbstractCompressibleMaterialLaw`. The `CompressibleMooneyRivlinMaterialLaw`
         * is one such example; instantiate one of these
         */
        CompressibleMooneyRivlinMaterialLaw<2> law(1.0, 0.5);

        /* For this problem, we fix the nodes on the surface for which Y < -0.9 */
        std::vector<unsigned> fixed_nodes;
        for ( TetrahedralMesh<2,2>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorBegin();
              iter != mesh.GetBoundaryNodeIteratorEnd();
              ++iter)
        {
            double Y = (*iter)->rGetLocation()[1];
            if(Y < -0.9)
            {
                fixed_nodes.push_back((*iter)->GetIndex());
            }
        }

        /* We will (later) apply Neumann boundary conditions to surface elements which lie below Y=0,
         * and these are collected here. (Minor, subtle, comment: we don't bother here checking Y>-0.9,
         * so the surface elements collected here include the ones on the Dirichlet boundary. This doesn't
         * matter as the Dirichlet boundary conditions to the nonlinear system and essentially overwrite
         * an Neumann related effects.
         */
        std::vector<BoundaryElement<1,2>*> boundary_elems;
        for (TetrahedralMesh<2,2>::BoundaryElementIterator iter
              = mesh.GetBoundaryElementIteratorBegin();
            iter != mesh.GetBoundaryElementIteratorEnd();
            ++iter)
        {
           BoundaryElement<1,2>* p_element = *iter;
           if(p_element->CalculateCentroid()[1]<0.0)
           {
               boundary_elems.push_back(p_element);
           }
        }
        assert(boundary_elems.size()>0);

        /* Create the problem definition class, and to start off with specify the zero
         * displacement nodes and an upward body force
         */
        SolidMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);

        c_vector<double,2> gravity;
        gravity(0) = 0;
        gravity(1) = 0.1;

        problem_defn.SetBodyForce(gravity);

        /* Declare the compressible solver, which has the same interface as the incompressible
         * one, and call `Solve()`
         */
        CompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                        problem_defn,
                                                        &law,
                                                        "CompressibleSolidMechanicsExample");
        solver.Solve();

        /* Now we call add additional boundary conditions, and call `Solve() again. Firstly: these
         * Neumnann conditions here are not specified traction boundary conditions (such BCs are specified
         * on the undeformed body), but instead, the (more natural) specification of a pressure
         * exactly in the ''normal direction on the deformed body''. We have to provide a set of boundary
         * elements of the mesh, and a pressure to act on those elements. The solver will automatically
         * compute the deformed normal directions on which the pressure acts. Note: with this type of
         * BC, the ordering of the nodes on the boundary elements needs to be consistent, otherwise some
         * normals will be computed to be inward and others outward. The solver will check this on the
         * original mesh and throw an exception if this is not the case. (The required ordering is such that:
         * in 2D, surface element nodes are ordered anticlockwise (looking at the whole mesh); in 3D, looking
         * at any surface element from outside the mesh, the three nodes are ordered anticlockwise. (Triangle/tetgen
         * automatically create boundary elements like this)).
         */
        double external_pressure = -0.4; // negative sign => inward pressure
        problem_defn.SetApplyNormalPressureOnDeformedSurface(boundary_elems, external_pressure);
        /* Call `Solve()` again, so now solving with fixed nodes, gravity, and pressure. The solution from
         * the previous solve will be used as the initial guess. Although at the moment the solution from the
         * previous call to `Solve()` will be over-written, calling `Solve()` repeatedly may be useful for
         * some problems: sometimes, Newton's method will fail to converge for given force/pressures etc, and can
         * be (very) helpful to ''increment'' the loading. For example, set the gravity to be (0,-9.81/3), solve,
         * then set it to be (0,-2*9.81/3), solve again, and finally set it to be (0,-9.81) and solve for a
         * final time
         */
        solver.Solve();
    }
};

#endif // TESTSOLVINGMOREELASTICITYPROBLEMSTUTORIALS_HPP_
