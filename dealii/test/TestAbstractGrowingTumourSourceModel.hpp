#ifndef TESTABSTRACTGROWINGTUMOURSOURCEMODEL_HPP_
#define TESTABSTRACTGROWINGTUMOURSOURCEMODEL_HPP_

#include <cxxtest/TestSuite.h>

#include "AbstractGrowingTumourSourceModel.hpp"
#include "SimpleTumourSourceModel.hpp"
#include "ConstantTumourSourceModel.hpp"
#include "FiniteElasticityAssembler.hpp"
#include "FiniteElasticityTools.hpp"
#include "MooneyRivlinMaterialLaw.hpp"
#include "TriangulationVertexIterator.hpp"

class TestAbstractGrowingTumourSourceModel : public CxxTest::TestSuite
{
public:
    void TestWithSimpleModel()
    {
        SimpleTumourSourceModel<2> simple_model;
        
        // add a couple of points to watch
        Point<2> position;
        simple_model.AddEvaluationPoint(0, position, 0);
        simple_model.AddEvaluationPoint(1, position, 0);
        
        TS_ASSERT_EQUALS(simple_model.GetNumEvaluationPoints(), 2);
        
        // the source value at the tumour points should be zero, 
        // as model hasn't been run yet
        TS_ASSERT_DELTA(simple_model.GetSourceValue(0), 0.0, 1e-12);
        TS_ASSERT_DELTA(simple_model.GetSourceValue(1), 0.0, 1e-12);        

        // run
        simple_model.Run(0,1);
        
        // this source model just return s = index,
        TS_ASSERT_DELTA(simple_model.GetSourceValue(0), 0.0, 1e-12);
        TS_ASSERT_DELTA(simple_model.GetSourceValue(1), 1.0, 1e-12);
        
        // exceptions:
        // evaluation point with index=0 already exists
        TS_ASSERT_THROWS_ANYTHING(simple_model.AddEvaluationPoint(0, position, 0));
        // no evaluation point corresponding to index=10 exists
        TS_ASSERT_THROWS_ANYTHING(simple_model.GetSourceValue(10));
    }

    void TestWithConstantModel()
    {
        double value = 0.343234;
        ConstantTumourSourceModel<2> constant_model(value);
        
        // add a points to watch
        Point<2> position;
        constant_model.AddEvaluationPoint(0, position, 0);

        // run
        constant_model.Run(0,1);
        
        // this source model just return s = index,
        TS_ASSERT_DELTA(constant_model.GetSourceValue(0), value, 1e-12);
    }

    void TestUpdateEvaluationPointsNewPosition()
    {
        Vector<double> body_force(2);
        body_force(0) = 1.0;
    
        MooneyRivlinMaterialLaw<2> mooney_rivlin_law(2.0);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0); 
        mesh.refine_global(1);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0);
        
        FiniteElasticityAssembler<2> finite_elasticity(&mesh,
                                                       &mooney_rivlin_law,
                                                       body_force,
                                                       1.0,
                                                       "temp");

        double value = 2.5;
        ConstantTumourSourceModel<2> constant_model(value);
        
        // add the first few nodes
        TriangulationVertexIterator<2> iter(&mesh); 

        //std::cout << "Adding " << iter.GetVertexGlobalIndex() << " - " << iter.GetVertex()[0] << " " <<iter.GetVertex()[1] << "\n";  
        constant_model.AddEvaluationPoint(0, iter.GetVertex(), iter.GetVertexGlobalIndex());
        iter.Next();

        //std::cout << "Adding " << iter.GetVertexGlobalIndex() << " - " << iter.GetVertex()[0] << " " <<iter.GetVertex()[1] << "\n";  
        constant_model.AddEvaluationPoint(1, iter.GetVertex(), iter.GetVertexGlobalIndex());

        //////////////////////////////////////////////////////////////////////////
        // check the old position and the new position of the eval points agree
        // Note: we are a friend class so can access the protect member variables
        //////////////////////////////////////////////////////////////////////////
        for(unsigned i=0; i<2; i++)
        {        
            Point<2> old_posn_pt = constant_model.mEvaluationPoints[i].OldPosition;
            Point<2> new_posn_pt = constant_model.mEvaluationPoints[i].NewPosition;
        
            TS_ASSERT_DELTA(old_posn_pt[0], new_posn_pt[0], 1e-12);
            TS_ASSERT_DELTA(old_posn_pt[1], new_posn_pt[1], 1e-12);
        }

        //////////////////////////////////////////////////////////////////////////
        // call update new position and check they still agree
        //////////////////////////////////////////////////////////////////////////
        constant_model.UpdateEvaluationPointsNewPosition(&finite_elasticity);
        for(unsigned i=0; i<2; i++)
        {        
            Point<2> old_posn_pt = constant_model.mEvaluationPoints[i].OldPosition;
            Point<2> new_posn_pt = constant_model.mEvaluationPoints[i].NewPosition;
        
            TS_ASSERT_DELTA(old_posn_pt[0], new_posn_pt[0], 1e-12);
            TS_ASSERT_DELTA(old_posn_pt[1], new_posn_pt[1], 1e-12);
        }
        
        //////////////////////////////////////////////////////////////////////////
        // solve the static system and update again
        //////////////////////////////////////////////////////////////////////////
        finite_elasticity.Solve();
        constant_model.UpdateEvaluationPointsNewPosition(&finite_elasticity);
        
        // evaluation point 0: equivalent to x=y=0, on fixed boundary, 
        // new point should be the same as the old point
        Point<2> old_posn_pt = constant_model.mEvaluationPoints[0].OldPosition;
        Point<2> new_posn_pt = constant_model.mEvaluationPoints[0].NewPosition;
        TS_ASSERT_DELTA(old_posn_pt[0], new_posn_pt[0], 1e-12);
        TS_ASSERT_DELTA(old_posn_pt[1], new_posn_pt[1], 1e-12);
          
        // evaluation point 1: check against hardcoded results
        old_posn_pt = constant_model.mEvaluationPoints[1].OldPosition;
        new_posn_pt = constant_model.mEvaluationPoints[1].NewPosition;
        TS_ASSERT_DELTA(old_posn_pt[0], 0.5,    1e-12);
        TS_ASSERT_DELTA(new_posn_pt[0], 0.5123, 1e-3);
        TS_ASSERT_DELTA(old_posn_pt[1], 0.0,    1e-12);
        TS_ASSERT_DELTA(new_posn_pt[1], 0.0144, 1e-3);

        ///////////////////////////////////////////
        // add another eval point and do it again
        ///////////////////////////////////////////       
        iter.Next();
        //std::cout << "Adding " << iter.GetVertexGlobalIndex() << " - " << iter.GetVertex()[0] << " " <<iter.GetVertex()[1] << "\n";  
        constant_model.AddEvaluationPoint(2, iter.GetVertex(), iter.GetVertexGlobalIndex());

        constant_model.UpdateEvaluationPointsNewPosition(&finite_elasticity);
        
        // evaluation point 0: equivalent to x=y=0, on fixed boundary, 
        // new point should be the same as the old point
        old_posn_pt = constant_model.mEvaluationPoints[0].OldPosition;
        new_posn_pt = constant_model.mEvaluationPoints[0].NewPosition;
        TS_ASSERT_DELTA(old_posn_pt[0], new_posn_pt[0], 1e-12);
        TS_ASSERT_DELTA(old_posn_pt[1], new_posn_pt[1], 1e-12);
          
        // evaluation point 1: check against hardcoded results
        old_posn_pt = constant_model.mEvaluationPoints[1].OldPosition;
        new_posn_pt = constant_model.mEvaluationPoints[1].NewPosition;
        TS_ASSERT_DELTA(old_posn_pt[0], 0.5,    1e-12);
        TS_ASSERT_DELTA(new_posn_pt[0], 0.5123, 1e-3);
        TS_ASSERT_DELTA(old_posn_pt[1], 0.0,    1e-12);
        TS_ASSERT_DELTA(new_posn_pt[1], 0.0144, 1e-3);
        
        old_posn_pt = constant_model.mEvaluationPoints[2].OldPosition;
        new_posn_pt = constant_model.mEvaluationPoints[2].NewPosition;
        TS_ASSERT_DELTA(old_posn_pt[0], 0.5,    1e-12);
        TS_ASSERT_DELTA(new_posn_pt[0], 0.5155, 1e-3);
        TS_ASSERT_DELTA(old_posn_pt[1], 0.5,    1e-12);
        TS_ASSERT_DELTA(new_posn_pt[1], 0.4999, 1e-3);
    } 
};


#endif /*TESTABSTRACTGROWINGTUMOURSOURCEMODEL_HPP_*/
