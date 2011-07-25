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


#ifndef TESTSOLIDMECHANICSPROBLEMDEFINITION_HPP_
#define TESTSOLIDMECHANICSPROBLEMDEFINITION_HPP_


#include <cxxtest/TestSuite.h>
#include "SolidMechanicsProblemDefinition.hpp"

c_vector<double,2> SomeFunction(c_vector<double,2>& rX, double t)
{
    c_vector<double,2> body_force;
    body_force(0) = rX(0)+t;
    body_force(1) = 2*(rX(1)+t);
    return body_force;
}

c_vector<double,2> AnotherFunction(c_vector<double,2>& rX, double t)
{
    c_vector<double,2> body_force;
    body_force(0) = rX(0)*t;
    body_force(1) = 10*rX(1)*t;
    return body_force;
}


class TestSolidMechanicsProblemDefinition : public CxxTest::TestSuite
{
public:
    void TestDefinition() throw(Exception)
    {
        QuadraticMesh<2> mesh(0.5, 1.0, 1.0);

        SolidMechanicsProblemDefinition<2> problem_defn(mesh);

        TS_ASSERT_DELTA(problem_defn.GetDensity(), 1.0, 1e-12);

        TS_ASSERT_EQUALS(problem_defn.GetBodyForceType(), CONSTANT_BODY_FORCE);
        TS_ASSERT_DELTA(problem_defn.GetConstantBodyForce()(0), 0.0, 1e-12);
        TS_ASSERT_DELTA(problem_defn.GetConstantBodyForce()(1), 0.0, 1e-12);

        TS_ASSERT_EQUALS(problem_defn.GetTractionBoundaryConditionType(), NO_TRACTIONS);

        problem_defn.SetDensity(2.0);
        TS_ASSERT_DELTA(problem_defn.GetDensity(), 2.0, 1e-12);


        //////////////////////////////////
        // Body force
        //////////////////////////////////

        problem_defn.SetBodyForce(SomeFunction);
        TS_ASSERT_EQUALS(problem_defn.GetBodyForceType(), FUNCTIONAL_BODY_FORCE);
        c_vector<double,2> X;
        X(0) = 10.0;
        X(1) = 11.0;
        double t = 0.5;
        TS_ASSERT_DELTA(problem_defn.EvaluateBodyForceFunction(X,t)(0), 10.5, 1e-12);
        TS_ASSERT_DELTA(problem_defn.EvaluateBodyForceFunction(X,t)(1), 23.0, 1e-12);

        c_vector<double,2> body_force;
        body_force(0) = -9.8;
        body_force(1) = 0.01;
        problem_defn.SetBodyForce(body_force);
        TS_ASSERT_EQUALS(problem_defn.GetBodyForceType(), CONSTANT_BODY_FORCE);
        TS_ASSERT_DELTA(problem_defn.GetConstantBodyForce()(0), -9.8,  1e-12);
        TS_ASSERT_DELTA(problem_defn.GetConstantBodyForce()(1),  0.01, 1e-12);


        //////////////////////////////////
        // Traction
        //////////////////////////////////

        std::vector<BoundaryElement<1,2>*> boundary_elements;
        std::vector<c_vector<double,2> > tractions;

        TetrahedralMesh<2,2>::BoundaryElementIterator iter
           = mesh.GetBoundaryElementIteratorBegin();

        c_vector<double,2> vec = zero_vector<double>(2);
        vec(0)=1.0;
        boundary_elements.push_back(*iter);
        tractions.push_back(vec);

        ++iter;
        vec(1)=2.0;
        boundary_elements.push_back(*iter);
        tractions.push_back(vec);

        problem_defn.SetTractionBoundaryConditions(boundary_elements, tractions);

        TS_ASSERT_EQUALS(problem_defn.GetTractionBoundaryConditionType(), ELEMENTWISE_TRACTION);

        TS_ASSERT_EQUALS(problem_defn.rGetTractionBoundaryElements().size(), 2u);
        TS_ASSERT_EQUALS(problem_defn.rGetElementwiseTractions().size(), 2u);

        // comparing addresses
        TS_ASSERT_EQUALS(problem_defn.rGetTractionBoundaryElements()[0], boundary_elements[0]);
        TS_ASSERT_EQUALS(problem_defn.rGetTractionBoundaryElements()[1], boundary_elements[1]);

        TS_ASSERT_DELTA(problem_defn.rGetElementwiseTractions()[0](0), 1.0, 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetElementwiseTractions()[0](1), 0.0, 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetElementwiseTractions()[1](0), 1.0, 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetElementwiseTractions()[1](1), 2.0, 1e-12);

        ++iter;
        boundary_elements.push_back(*iter);
        double pressure = 3423.342;

        problem_defn.SetApplyNormalPressureOnDeformedSurface(boundary_elements, pressure);

        TS_ASSERT_EQUALS(problem_defn.rGetTractionBoundaryElements().size(), 3u);

        // comparing addresses
        TS_ASSERT_EQUALS(problem_defn.rGetTractionBoundaryElements()[0], boundary_elements[0]);
        TS_ASSERT_EQUALS(problem_defn.rGetTractionBoundaryElements()[1], boundary_elements[1]);
        TS_ASSERT_EQUALS(problem_defn.rGetTractionBoundaryElements()[2], boundary_elements[2]);

        TS_ASSERT_DELTA(problem_defn.GetNormalPressure(), 3423.342, 1e-12);

        ++iter;
        boundary_elements.push_back(*iter);
        problem_defn.SetTractionBoundaryConditions(boundary_elements, AnotherFunction);
        TS_ASSERT_EQUALS(problem_defn.rGetTractionBoundaryElements().size(), 4u);

        TS_ASSERT_DELTA(problem_defn.EvaluateTractionFunction(X,t)(0), 5.0,  1e-12);
        TS_ASSERT_DELTA(problem_defn.EvaluateTractionFunction(X,t)(1), 55.0, 1e-12);


        //////////////////////////////////
        // Fixed nodes
        //////////////////////////////////

        std::vector<unsigned> fixed_nodes;
        fixed_nodes.push_back(0);
        fixed_nodes.push_back(4);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        TS_ASSERT_EQUALS(problem_defn.rGetFixedNodes().size(), 2u);
        TS_ASSERT_EQUALS(problem_defn.rGetFixedNodes()[0], 0u);
        TS_ASSERT_EQUALS(problem_defn.rGetFixedNodes()[1], 4u);
        TS_ASSERT_EQUALS(problem_defn.rGetFixedNodeDisplacements().size(), 2u);

        TS_ASSERT_DELTA(problem_defn.rGetFixedNodeDisplacements()[0](0), 0.0, 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetFixedNodeDisplacements()[0](1), 0.0, 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetFixedNodeDisplacements()[1](0), 0.0, 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetFixedNodeDisplacements()[1](1), 0.0, 1e-12);

        fixed_nodes.push_back(8);

        std::vector<c_vector<double,2> > locations;
        c_vector<double,2> location = zero_vector<double>(2);
        locations.push_back(location);
        location(1)=0.1;
        locations.push_back(location);
        location(0)=0.1;
        locations.push_back(location);

        problem_defn.SetFixedNodes(fixed_nodes, locations);

        TS_ASSERT_EQUALS(problem_defn.rGetFixedNodes().size(), 3u);
        TS_ASSERT_EQUALS(problem_defn.rGetFixedNodes()[0], 0u);
        TS_ASSERT_EQUALS(problem_defn.rGetFixedNodes()[1], 4u);
        TS_ASSERT_EQUALS(problem_defn.rGetFixedNodes()[2], 8u);
        TS_ASSERT_EQUALS(problem_defn.rGetFixedNodeDisplacements().size(), 3u);

        TS_ASSERT_DELTA(problem_defn.rGetFixedNodeDisplacements()[0](0), 0.0 - mesh.GetNode(0)->rGetLocation()[0], 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetFixedNodeDisplacements()[0](1), 0.0 - mesh.GetNode(0)->rGetLocation()[1], 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetFixedNodeDisplacements()[1](0), 0.0 - mesh.GetNode(4)->rGetLocation()[0], 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetFixedNodeDisplacements()[1](1), 0.1 - mesh.GetNode(4)->rGetLocation()[1], 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetFixedNodeDisplacements()[2](0), 0.1 - mesh.GetNode(8)->rGetLocation()[0], 1e-12);
        TS_ASSERT_DELTA(problem_defn.rGetFixedNodeDisplacements()[2](1), 0.1 - mesh.GetNode(8)->rGetLocation()[1], 1e-12);

    }
};

#endif /* TESTSOLIDMECHANICSPROBLEMDEFINITION_HPP_ */
