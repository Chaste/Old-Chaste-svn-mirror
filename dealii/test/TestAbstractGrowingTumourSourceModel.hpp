#ifndef TESTABSTRACTGROWINGTUMOURSOURCEMODEL_HPP_
#define TESTABSTRACTGROWINGTUMOURSOURCEMODEL_HPP_

#include <cxxtest/TestSuite.h>

#include "AbstractGrowingTumourSourceModel.hpp"
#include "SimpleTumourSourceModel.hpp"
#include "ConstantTumourSourceModel.hpp"

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


};


#endif /*TESTABSTRACTGROWINGTUMOURSOURCEMODEL_HPP_*/
