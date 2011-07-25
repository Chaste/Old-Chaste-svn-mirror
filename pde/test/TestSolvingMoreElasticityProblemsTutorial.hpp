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

/* Comments to be added .. */

#include <cxxtest/TestSuite.h>
#include "CompressibleNonlinearElasticitySolver.hpp"
#include "CompressibleMooneyRivlinMaterialLaw.hpp"
#include "TrianglesMeshReader.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestSolvingMoreElasticityProblemsTutorials : public CxxTest::TestSuite
{
public:
    void TestSolvingCompressibleProblem() throw (Exception)
    {
        QuadraticMesh<2> mesh;
        TrianglesMeshReader<2,2> reader("mesh/test/data/disk_522_elements");
        mesh.ConstructFromLinearMeshReader(reader);

        CompressibleMooneyRivlinMaterialLaw<2> law(1.0, 0.5);

        /* For fixed nodes, we choose the nodes for which Y
         *
         */
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

        SolidMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);

        c_vector<double,2> gravity;
        gravity(0) = 0;
        gravity(1) = 0.1;

        problem_defn.SetBodyForce(gravity);

        CompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                        problem_defn,
                                                        &law,
                                                        "CompressibleSolidMechanicsExample");
        solver.Solve();


        double external_pressure = -0.4;
        problem_defn.SetApplyNormalPressureOnDeformedSurface(boundary_elems, external_pressure);

        solver.Solve();
    }
};

#endif // TESTSOLVINGMOREELASTICITYPROBLEMSTUTORIALS_HPP_
