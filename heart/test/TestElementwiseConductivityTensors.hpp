#ifndef TESTELEMENTWISECONDUCTIVITYTENSORS_HPP_
#define TESTELEMENTWISECONDUCTIVITYTENSORS_HPP_

#include <cxxtest/TestSuite.h>
#include "ElementwiseConductivityTensors.hpp"

class TestElementwiseConductivityTensors : public CxxTest::TestSuite
{
public:
    void TestConstantTensor3D()
    {
        double longitudinal_const = 2.1, transverse_const = 0.8, normal_const = 0.135;        
        
        ElementwiseConductivityTensors<3> simple_tensors;
        simple_tensors.SetConstantConductivities(longitudinal_const, transverse_const, normal_const);
        simple_tensors.Init();
        
        for (unsigned tensor_index=0; tensor_index<4; tensor_index++)
        {
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](0,0), longitudinal_const);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](0,1), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](0,2), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](1,0), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](1,1), transverse_const);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](1,2), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](2,0), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](2,1), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](2,2), normal_const);
        }    
    }

    void TestConstantTensor1D()
    {
        double longitudinal_const = 2.1;        
        
        ElementwiseConductivityTensors<1> simple_tensors;
        simple_tensors.SetConstantConductivities(longitudinal_const);
        simple_tensors.Init();
        
        for (unsigned tensor_index=0; tensor_index<4; tensor_index++)
        {
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](0,0), longitudinal_const);
        }    
    }


    void TestFibreOrientationTensor3D()
    {   
    	double longitudinal_const = 2.1, transverse_const = 0.8, normal_const = 0.135;        
        
        ElementwiseConductivityTensors<3> simple_tensors;
        simple_tensors.SetConstantConductivities(longitudinal_const, transverse_const, normal_const);
    	simple_tensors.SetFibreOrientationFile("heart/test/data/SimpleFibreOrientation3D.fibres");
        simple_tensors.Init();
        
        for (unsigned tensor_index=0; tensor_index<4; tensor_index++)
        {
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](0,0), longitudinal_const);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](0,1), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](0,2), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](1,0), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](1,1), transverse_const);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](1,2), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](2,0), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](2,1), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](2,2), normal_const);
        }    
    }
    
    void TestHeterogeneousConductivitiesTensor3D()
    {
        std::vector<double> longitudinal_cond, transverse_cond, normal_cond;
        
        longitudinal_cond.push_back(0.0);
        longitudinal_cond.push_back(100.0);
        longitudinal_cond.push_back(200.0);
        longitudinal_cond.push_back(300.0);
        
        transverse_cond.push_back(0.0);
        transverse_cond.push_back(10.0);
        transverse_cond.push_back(20.0);
        transverse_cond.push_back(30.0);
        
        normal_cond.push_back(0.0);
        normal_cond.push_back(1.0);
        normal_cond.push_back(2.0);
        normal_cond.push_back(3.0);
        
        ElementwiseConductivityTensors<3> simple_tensors;
        simple_tensors.SetNonConstantConductivities(&longitudinal_cond, &transverse_cond, &normal_cond);
        simple_tensors.Init();
        
        for (unsigned tensor_index=0; tensor_index<4; tensor_index++)
        {
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](0,0), 100*tensor_index);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](0,1), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](0,2), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](1,0), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](1,1), 10*tensor_index);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](1,2), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](2,0), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](2,1), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](2,2), tensor_index);
        }    
    }

    void TestHeterogeneousCondPlusFibreOrientationTensor3D()
    {
        std::vector<double> longitudinal_cond;
        std::vector<double> transverse_cond;
        std::vector<double> normal_cond;
        
        longitudinal_cond.push_back(0.0);
        longitudinal_cond.push_back(100.0);
        longitudinal_cond.push_back(200.0);
        longitudinal_cond.push_back(300.0);
        
        transverse_cond.push_back(0.0);
        transverse_cond.push_back(10.0);
        transverse_cond.push_back(20.0);
        transverse_cond.push_back(30.0);
        
        normal_cond.push_back(0.0);
        normal_cond.push_back(1.0);
        normal_cond.push_back(2.0);
        normal_cond.push_back(3.0);
        
        ElementwiseConductivityTensors<3> simple_tensors;
        simple_tensors.SetNonConstantConductivities(&longitudinal_cond, &transverse_cond, &normal_cond);
        simple_tensors.SetFibreOrientationFile("heart/test/data/SimpleFibreOrientation3D.fibres");
        simple_tensors.Init();
        
        for (unsigned tensor_index=0; tensor_index<4; tensor_index++)
        {
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](0,0), 100*tensor_index);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](0,1), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](0,2), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](1,0), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](1,1), 10*tensor_index);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](1,2), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](2,0), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](2,1), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](2,2), tensor_index);
        }    
    }

    void TestHeterogeneousCondPlusFibreOrientationTensor2D()
    {
        std::vector<double> longitudinal_cond;
        std::vector<double> transverse_cond;
        
        longitudinal_cond.push_back(0.0);
        longitudinal_cond.push_back(100.0);
        longitudinal_cond.push_back(200.0);
        longitudinal_cond.push_back(300.0);
        
        transverse_cond.push_back(0.0);
        transverse_cond.push_back(10.0);
        transverse_cond.push_back(20.0);
        transverse_cond.push_back(30.0);
        
        ElementwiseConductivityTensors<2> simple_tensors;
        simple_tensors.SetNonConstantConductivities(&longitudinal_cond, &transverse_cond);
        simple_tensors.SetFibreOrientationFile("heart/test/data/SimpleFibreOrientation2D.fibres");
        simple_tensors.Init();
        
        for (unsigned tensor_index=0; tensor_index<4; tensor_index++)
        {
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](0,0), 100*tensor_index);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](0,1), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](1,0), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](1,1), 10*tensor_index);
        }    
    }

    void TestHeterogeneousCondPlusFibreOrientationTensor1D()
    {
        std::vector<double> longitudinal_cond;
        
        longitudinal_cond.push_back(0.0);
        longitudinal_cond.push_back(100.0);
        longitudinal_cond.push_back(200.0);
        longitudinal_cond.push_back(300.0);
        
        ElementwiseConductivityTensors<1> simple_tensors;
        simple_tensors.SetNonConstantConductivities(&longitudinal_cond);
        simple_tensors.SetFibreOrientationFile("heart/test/data/SimpleFibreOrientation1D.fibres");
        simple_tensors.Init();
        
        for (unsigned tensor_index=0; tensor_index<4; tensor_index++)
        {
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](0,0), 100*tensor_index);
        }    
    }

};

#endif /*TESTFIBREORIENTATIONTENSORS_HPP_*/
