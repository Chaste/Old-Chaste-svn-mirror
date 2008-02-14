#ifndef TESTELEMENTWISECONDUCTIVITYTENSORS_HPP_
#define TESTELEMENTWISECONDUCTIVITYTENSORS_HPP_

#include "UblasCustomFunctions.hpp"

#include <cxxtest/TestSuite.h>
#include "ElementwiseConductivityTensors.hpp"


class TestElementwiseConductivityTensors : public CxxTest::TestSuite
{
public:
    void TestConstantTensor3D()
    {
        ElementwiseConductivityTensors<3> simple_tensors;               
        simple_tensors.SetConstantConductivities(Create_c_vector(2.1, 0.8, 0.135));
        simple_tensors.Init();
        
        for (unsigned tensor_index=0; tensor_index<4; tensor_index++)
        {
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](0,0), 2.1);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](0,1), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](0,2), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](1,0), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](1,1), 0.8);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](1,2), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](2,0), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](2,1), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](2,2), 0.135);
        }    
    }

    void TestTensorException() throw (Exception)
    {
        c_vector<double, 2> constant_conductivities(Create_c_vector(2.1, 0.8));                       
        
        ElementwiseConductivityTensors<1> simple_1d_tensors;
        TS_ASSERT_THROWS_ANYTHING(simple_1d_tensors.SetConstantConductivities(constant_conductivities));                     

        ElementwiseConductivityTensors<3> simple_3d_tensors;
        TS_ASSERT_THROWS_ANYTHING(simple_3d_tensors.SetConstantConductivities(constant_conductivities));                     
    }

    void TestFibreOrientationTensor3D()
    {   
        c_vector<double, 3> constant_conductivities(Create_c_vector(2.1,0.8,0.135));
        
        ElementwiseConductivityTensors<3> simple_tensors;
        simple_tensors.SetConstantConductivities(constant_conductivities);
    	simple_tensors.SetFibreOrientationFile("heart/test/data/SimpleFibreOrientation3D.fibres");
        simple_tensors.Init();
        
        for (unsigned tensor_index=0; tensor_index<4; tensor_index++)
        {
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](0,0), constant_conductivities[0]);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](0,1), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](0,2), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](1,0), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](1,1), constant_conductivities[1]);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](1,2), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](2,0), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](2,1), 0.0);
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](2,2), constant_conductivities[2]);
        }    
    }
    
    void TestHeterogeneousConductivitiesTensor3D()
    {
        std::vector<c_vector<double, 3> > non_constant_conductivities;
        non_constant_conductivities.push_back(Create_c_vector(0,0,0));
        non_constant_conductivities.push_back(Create_c_vector(100,10,1));
        non_constant_conductivities.push_back(Create_c_vector(200,20,2));
        non_constant_conductivities.push_back(Create_c_vector(300,30,3));
        
        ElementwiseConductivityTensors<3> simple_tensors;
        simple_tensors.SetNonConstantConductivities(&non_constant_conductivities);
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

    void TestHeterogeneousCondPlusFibreOrientationTensor1D()
    {
        std::vector<c_vector<double, 1> > non_constant_conductivities;
        non_constant_conductivities.push_back(Create_c_vector(0));
        non_constant_conductivities.push_back(Create_c_vector(100));
        non_constant_conductivities.push_back(Create_c_vector(200));
        non_constant_conductivities.push_back(Create_c_vector(300));
        
        ElementwiseConductivityTensors<1> simple_tensors;
        simple_tensors.SetNonConstantConductivities(&non_constant_conductivities);
        simple_tensors.SetFibreOrientationFile("heart/test/data/SimpleFibreOrientation1D.fibres");
        simple_tensors.Init();
        
        for (unsigned tensor_index=0; tensor_index<4; tensor_index++)
        {
            TS_ASSERT_EQUALS(simple_tensors[tensor_index](0,0), 100*tensor_index);
        }    
    }

    void TestHeterogeneousCondPlusFibreOrientationTensor2D()
    {
        std::vector<c_vector<double, 2> > non_constant_conductivities;
        non_constant_conductivities.push_back(Create_c_vector(0,0));
        non_constant_conductivities.push_back(Create_c_vector(100,10));
        non_constant_conductivities.push_back(Create_c_vector(200,20));
        non_constant_conductivities.push_back(Create_c_vector(300,30));
        
        ElementwiseConductivityTensors<2> simple_tensors;
        simple_tensors.SetNonConstantConductivities(&non_constant_conductivities);
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

    void TestHeterogeneousCondPlusFibreOrientationTensor3D()
    {
        std::vector<c_vector<double, 3> > non_constant_conductivities;
        non_constant_conductivities.push_back(Create_c_vector(0,0,0));
        non_constant_conductivities.push_back(Create_c_vector(100,10,1));
        non_constant_conductivities.push_back(Create_c_vector(200,20,2));
        non_constant_conductivities.push_back(Create_c_vector(300,30,3));
        
        ElementwiseConductivityTensors<3> simple_tensors;
        simple_tensors.SetNonConstantConductivities(&non_constant_conductivities);
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

};

#endif /*TESTFIBREORIENTATIONTENSORS_HPP_*/
