/*

Copyright (C) University of Oxford, 2005-2009

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


#ifndef TESTABSTRACTGROWINGTUMOURSOURCEMODEL_HPP_
#define TESTABSTRACTGROWINGTUMOURSOURCEMODEL_HPP_

#include <cxxtest/TestSuite.h>

#include "AbstractGrowingTumourSourceModel.hpp"
#include "SimpleTumourSourceModel.hpp"
#include "ConstantTumourSourceModel.hpp"
#include "FiniteElasticityAssembler.hpp"
#include "FiniteElasticityTools.hpp"
#include "MooneyRivlinMaterialLaw2.hpp"
#include "TriangulationVertexIterator.hpp"

class TestAbstractGrowingTumourSourceModel : public CxxTest::TestSuite
{
public:
    void TestWithSimpleModel()
    {
        SimpleTumourSourceModel<2> simple_model;

        // add a couple of points to watch
        Point<2> position;
        simple_model.AddEvaluationPoint(0, position);
        simple_model.AddEvaluationPoint(1, position);

        TS_ASSERT_EQUALS(simple_model.GetNumEvaluationPoints(), 2u);

        // the source value at the tumour points should be zero,
        // as model hasn't been run yet
        TS_ASSERT_DELTA(simple_model.GetSourceValue(0), 0.0, 1e-12);
        TS_ASSERT_DELTA(simple_model.GetSourceValue(1), 0.0, 1e-12);

        // run
        simple_model.Run(0,1,NULL);

        // this source model just return s = index,
        TS_ASSERT_DELTA(simple_model.GetSourceValue(0), 0.0, 1e-12);
        TS_ASSERT_DELTA(simple_model.GetSourceValue(1), 1.0, 1e-12);

        // exceptions:
        // evaluation point with index=0 already exists
        TS_ASSERT_THROWS_ANYTHING(simple_model.AddEvaluationPoint(0, position));
        // no evaluation point corresponding to index=10 exists
        TS_ASSERT_THROWS_ANYTHING(simple_model.GetSourceValue(10));
    }

    void TestWithConstantModel()
    {
        double value = 0.343234;
        ConstantTumourSourceModel<2> constant_model(value);

        // add a points to watch
        Point<2> position;
        constant_model.AddEvaluationPoint(0, position);

        // run
        constant_model.Run(0,1,NULL);

        // this source model just return s = index,
        TS_ASSERT_DELTA(constant_model.GetSourceValue(0), value, 1e-12);
    }

    void TestUpdateEvaluationPointsNewPosition()
    {
        Vector<double> body_force(2);
        body_force(0) = 1.0;

        MooneyRivlinMaterialLaw2<2> mooney_rivlin_law(2.0);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(1);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 0.0);

        FiniteElasticityAssembler<2> finite_elasticity(&mesh,
                                                       &mooney_rivlin_law,
                                                       body_force,
                                                       1.0,
                                                       "");

        double value = 2.5;
        ConstantTumourSourceModel<2> constant_model(value);

        // add the first few nodes

        TriangulationVertexIterator<2> iter(&mesh);

        //std::cout << "Adding " << iter.GetVertexGlobalIndex() << " - " << iter.GetVertex()[0] << " " <<iter.GetVertex()[1] << "\n";
        unsigned indices[3];
        indices[0] = iter.GetVertexGlobalIndex();
        constant_model.AddEvaluationPoint(indices[0], iter.GetVertex());

        iter.Next();

        //std::cout << "Adding " << iter.GetVertexGlobalIndex() << " - " << iter.GetVertex()[0] << " " <<iter.GetVertex()[1] << "\n";
        indices[1] = iter.GetVertexGlobalIndex();
        constant_model.AddEvaluationPoint(indices[1], iter.GetVertex());

        //////////////////////////////////////////////////////////////////////////
        // check the old position and the new position of the eval points agree
        // Note: we are a friend class so can access the protect member variables
        //////////////////////////////////////////////////////////////////////////
        for (unsigned i=0; i<2; i++)
        {
            Point<2> old_posn_pt = constant_model.mEvaluationPoints[indices[i]].OldPosition;
            Point<2> new_posn_pt = constant_model.mEvaluationPoints[indices[i]].NewPosition;

            TS_ASSERT_DELTA(old_posn_pt[0], new_posn_pt[0], 1e-12);
            TS_ASSERT_DELTA(old_posn_pt[1], new_posn_pt[1], 1e-12);
        }

        //////////////////////////////////////////////////////////////////////////
        // call update new position and check they still agree
        //////////////////////////////////////////////////////////////////////////
        constant_model.UpdateEvaluationPointsNewPosition(finite_elasticity.rGetDeformedPosition());
        for (unsigned i=0; i<2; i++)
        {
            Point<2> old_posn_pt = constant_model.mEvaluationPoints[indices[i]].OldPosition;
            Point<2> new_posn_pt = constant_model.mEvaluationPoints[indices[i]].NewPosition;

            TS_ASSERT_DELTA(old_posn_pt[0], new_posn_pt[0], 1e-12);
            TS_ASSERT_DELTA(old_posn_pt[1], new_posn_pt[1], 1e-12);
        }

        //////////////////////////////////////////////////////////////////////////
        // solve the static system and update again
        //////////////////////////////////////////////////////////////////////////
        finite_elasticity.StaticSolve();
        constant_model.UpdateEvaluationPointsNewPosition(finite_elasticity.rGetDeformedPosition());

        // evaluation point 0: equivalent to x=y=0, on fixed boundary,
        // new point should be the same as the old point
        Point<2> old_posn_pt = constant_model.mEvaluationPoints[indices[0]].OldPosition;
        Point<2> new_posn_pt = constant_model.mEvaluationPoints[indices[0]].NewPosition;
        TS_ASSERT_DELTA(old_posn_pt[0], new_posn_pt[0], 1e-12);
        TS_ASSERT_DELTA(old_posn_pt[1], new_posn_pt[1], 1e-12);

        // evaluation point 1: check against hardcoded results
        old_posn_pt = constant_model.mEvaluationPoints[indices[1]].OldPosition;
        new_posn_pt = constant_model.mEvaluationPoints[indices[1]].NewPosition;
        TS_ASSERT_DELTA(old_posn_pt[0], 0.5,    1e-12);
        TS_ASSERT_DELTA(new_posn_pt[0], 0.5123, 1e-3);
        TS_ASSERT_DELTA(old_posn_pt[1], 0.0,    1e-12);
        TS_ASSERT_DELTA(new_posn_pt[1], 0.0144, 1e-3);

        ///////////////////////////////////////////
        // add another eval point and do it again
        ///////////////////////////////////////////
        iter.Next();

        indices[2] = iter.GetVertexGlobalIndex();
        //std::cout << "Adding " << iter.GetVertexGlobalIndex() << " - " << iter.GetVertex()[0] << " " <<iter.GetVertex()[1] << "\n";

        constant_model.AddEvaluationPoint(indices[2], iter.GetVertex());

        constant_model.UpdateEvaluationPointsNewPosition(finite_elasticity.rGetDeformedPosition());

        // evaluation point 0: equivalent to x=y=0, on fixed boundary,
        // new point should be the same as the old point
        old_posn_pt = constant_model.mEvaluationPoints[indices[0]].OldPosition;
        new_posn_pt = constant_model.mEvaluationPoints[indices[0]].NewPosition;
        TS_ASSERT_DELTA(old_posn_pt[0], new_posn_pt[0], 1e-12);
        TS_ASSERT_DELTA(old_posn_pt[1], new_posn_pt[1], 1e-12);

        // evaluation point 1: check against hardcoded results
        old_posn_pt = constant_model.mEvaluationPoints[indices[1]].OldPosition;
        new_posn_pt = constant_model.mEvaluationPoints[indices[1]].NewPosition;
        TS_ASSERT_DELTA(old_posn_pt[0], 0.5,    1e-12);
        TS_ASSERT_DELTA(new_posn_pt[0], 0.5123, 1e-3);
        TS_ASSERT_DELTA(old_posn_pt[1], 0.0,    1e-12);
        TS_ASSERT_DELTA(new_posn_pt[1], 0.0144, 1e-3);

        old_posn_pt = constant_model.mEvaluationPoints[indices[2]].OldPosition;
        new_posn_pt = constant_model.mEvaluationPoints[indices[2]].NewPosition;
        TS_ASSERT_DELTA(old_posn_pt[0], 0.5,    1e-12);
        TS_ASSERT_DELTA(new_posn_pt[0], 0.5155, 1e-3);
        TS_ASSERT_DELTA(old_posn_pt[1], 0.5,    1e-12);
        TS_ASSERT_DELTA(new_posn_pt[1], 0.4999, 1e-3);
    }
};


#endif /*TESTABSTRACTGROWINGTUMOURSOURCEMODEL_HPP_*/
