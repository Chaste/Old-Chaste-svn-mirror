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
#ifndef TESTCONDUCTIVITYTENSORS_HPP_
#define TESTCONDUCTIVITYTENSORS_HPP_

#include "UblasCustomFunctions.hpp"

#include <cxxtest/TestSuite.h>
#include "OrthotropicConductivityTensors.hpp"
#include "AxisymmetricConductivityTensors.hpp"


class TestConductivityTensors : public CxxTest::TestSuite
{
public:
    void TestConstantTensor3D()
    {
        OrthotropicConductivityTensors<3> ortho_tensors;
        ortho_tensors.SetConstantConductivities(Create_c_vector(2.1, 0.8, 0.135));
        ortho_tensors.Init();

        for (unsigned tensor_index=0; tensor_index<4; tensor_index++)
        {
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](0,0), 2.1);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](0,1), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](0,2), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](1,0), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](1,1), 0.8);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](1,2), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](2,0), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](2,1), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](2,2), 0.135);
        }

        AxisymmetricConductivityTensors<3> axi_tensors;
        axi_tensors.SetConstantConductivities(Create_c_vector(2.1, 0.8, 0.8));
        axi_tensors.Init();

        for (unsigned tensor_index=0; tensor_index<4; tensor_index++)
        {
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](0,0), 2.1);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](0,1), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](0,2), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](1,0), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](1,1), 0.8);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](1,2), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](2,0), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](2,1), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](2,2), 0.8);
        }
    }

    void TestTensorException() throw (Exception)
    {
        c_vector<double, 2> constant_conductivities(Create_c_vector(2.1, 0.8));

        OrthotropicConductivityTensors<1> ortho_1d_tensors;
        TS_ASSERT_THROWS_ANYTHING(ortho_1d_tensors.SetConstantConductivities(constant_conductivities));

        OrthotropicConductivityTensors<3> ortho_3d_tensors;
        TS_ASSERT_THROWS_ANYTHING(ortho_3d_tensors.SetConstantConductivities(constant_conductivities));

        // AxisymmetricConductivityTensors only makes sense in 3D problems
        TS_ASSERT_THROWS_ANYTHING(AxisymmetricConductivityTensors<2> axi_tensor);

        // Transversal and longitudinal conductivities should have the same value
        AxisymmetricConductivityTensors<3> axi_3d_tensor;
        TS_ASSERT_THROWS_ANYTHING( axi_3d_tensor.SetConstantConductivities(Create_c_vector(0.5,0.25,0.15)) );
    }

    void TestFibreOrientationFileExceptions() throw (Exception)
    {
        OrthotropicConductivityTensors<3> ortho_tensors;
        ortho_tensors.SetConstantConductivities( Create_c_vector(2.1, 0.8, 0.135) );
        ortho_tensors.SetFibreOrientationFile("non_existing_file.fibres");
        TS_ASSERT_THROWS_ANYTHING(ortho_tensors.Init()); // non existing file

        ortho_tensors.SetFibreOrientationFile ("heart/test/data/SimpleOrthotropic2D.fibres");
        TS_ASSERT_THROWS_ANYTHING(ortho_tensors.Init()); // mismatching SPACE_DIM and # vectors in file

        ortho_tensors.SetFibreOrientationFile("heart/test/data/SimpleOrthotropic2DWrongFormat.fibres");
        TS_ASSERT_THROWS_ANYTHING(ortho_tensors.Init()); // wrong file format

        ortho_tensors.SetFibreOrientationFile("heart/test/data/SimpleOrthotropic3DShortFile.fibres");
        TS_ASSERT_THROWS_ANYTHING(ortho_tensors.Init()); // short file
    }

    void TestFibreOrientationTensor3D()
    {
        c_vector<double, 3> constant_conductivities(Create_c_vector(2.1,0.8,0.135));

        OrthotropicConductivityTensors<3> ortho_tensors;
        ortho_tensors.SetConstantConductivities(constant_conductivities);
        ortho_tensors.SetFibreOrientationFile("heart/test/data/SimpleOrthotropic3D.fibres");
        ortho_tensors.Init();

        for (unsigned tensor_index=0; tensor_index<4; tensor_index++)
        {
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](0,0), constant_conductivities[0]);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](0,1), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](0,2), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](1,0), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](1,1), constant_conductivities[1]);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](1,2), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](2,0), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](2,1), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](2,2), constant_conductivities[2]);
        }
    }

    void TestFibreOrientationAxisymmetric3D()
    {
        c_vector<double, 3> constant_conductivities(Create_c_vector(2.1,0.8,0.8));

        OrthotropicConductivityTensors<3> ortho_tensors;
        ortho_tensors.SetConstantConductivities(constant_conductivities);
        ortho_tensors.SetFibreOrientationFile("heart/test/data/SimpleAxisymmetric.fibre");
        TS_ASSERT_THROWS_ANYTHING(ortho_tensors.Init());

        AxisymmetricConductivityTensors<3> axi_tensors;
        axi_tensors.SetConstantConductivities(constant_conductivities);
        axi_tensors.SetFibreOrientationFile("heart/test/data/SimpleAxisymmetric.fibres");
        axi_tensors.Init();

        for (unsigned tensor_index=0; tensor_index<4; tensor_index++)
        {
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](0,0), constant_conductivities[0]);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](0,1), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](0,2), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](1,0), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](1,1), constant_conductivities[1]);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](1,2), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](2,0), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](2,1), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](2,2), constant_conductivities[1]);
        }
    }

    void TestHeterogeneousConductivitiesTensor3D()
    {
        std::vector<c_vector<double, 3> > non_constant_conductivities;
        non_constant_conductivities.push_back(Create_c_vector(0,0,0));
        non_constant_conductivities.push_back(Create_c_vector(100,10,1));
        non_constant_conductivities.push_back(Create_c_vector(200,20,2));
        non_constant_conductivities.push_back(Create_c_vector(300,30,3));

        OrthotropicConductivityTensors<3> ortho_tensors;
        ortho_tensors.SetNonConstantConductivities(&non_constant_conductivities);
        ortho_tensors.Init();

        for (unsigned tensor_index=0; tensor_index<4; tensor_index++)
        {
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](0,0), 100*tensor_index);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](0,1), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](0,2), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](1,0), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](1,1), 10*tensor_index);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](1,2), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](2,0), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](2,1), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](2,2), tensor_index);
        }

        AxisymmetricConductivityTensors<3> axi_tensors;
        axi_tensors.SetNonConstantConductivities(&non_constant_conductivities);
        axi_tensors.Init();

        for (unsigned tensor_index=0; tensor_index<4; tensor_index++)
        {
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](0,0), 100*tensor_index);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](0,1), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](0,2), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](1,0), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](1,1), 10*tensor_index);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](1,2), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](2,0), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](2,1), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](2,2), 10*tensor_index);
        }


    }

    void TestHeterogeneousCondPlusFibreOrientationTensor1D()
    {
        std::vector<c_vector<double, 1> > non_constant_conductivities;
        non_constant_conductivities.push_back(Create_c_vector(0));
        non_constant_conductivities.push_back(Create_c_vector(100));
        non_constant_conductivities.push_back(Create_c_vector(200));
        non_constant_conductivities.push_back(Create_c_vector(300));

        OrthotropicConductivityTensors<1> ortho_tensors;
        ortho_tensors.SetNonConstantConductivities(&non_constant_conductivities);
        ortho_tensors.SetFibreOrientationFile("heart/test/data/SimpleOrthotropic1D.fibres");
        ortho_tensors.Init();

        for (unsigned tensor_index=0; tensor_index<4; tensor_index++)
        {
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](0,0), 100*tensor_index);
        }
    }

    void TestHeterogeneousCondPlusFibreOrientationTensor2D()
    {
        std::vector<c_vector<double, 2> > non_constant_conductivities;
        non_constant_conductivities.push_back(Create_c_vector(0,0));
        non_constant_conductivities.push_back(Create_c_vector(100,10));
        non_constant_conductivities.push_back(Create_c_vector(200,20));
        non_constant_conductivities.push_back(Create_c_vector(300,30));

        OrthotropicConductivityTensors<2> ortho_tensors;
        ortho_tensors.SetNonConstantConductivities(&non_constant_conductivities);
        ortho_tensors.SetFibreOrientationFile("heart/test/data/SimpleOrthotropic2D.fibres");
        ortho_tensors.Init();

        for (unsigned tensor_index=0; tensor_index<4; tensor_index++)
        {
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](0,0), 100*tensor_index);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](0,1), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](1,0), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](1,1), 10*tensor_index);
        }
    }

    void TestHeterogeneousCondPlusFibreOrientationTensor3D()
    {
        std::vector<c_vector<double, 3> > non_constant_conductivities;
        non_constant_conductivities.push_back(Create_c_vector(0,0,0));
        non_constant_conductivities.push_back(Create_c_vector(100,10,1));
        non_constant_conductivities.push_back(Create_c_vector(200,20,2));
        non_constant_conductivities.push_back(Create_c_vector(300,30,3));

        OrthotropicConductivityTensors<3> ortho_tensors;
        ortho_tensors.SetNonConstantConductivities(&non_constant_conductivities);
        ortho_tensors.SetFibreOrientationFile("heart/test/data/SimpleOrthotropic3D.fibres");
        ortho_tensors.Init();

        for (unsigned tensor_index=0; tensor_index<4; tensor_index++)
        {
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](0,0), 100*tensor_index);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](0,1), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](0,2), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](1,0), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](1,1), 10*tensor_index);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](1,2), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](2,0), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](2,1), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](2,2), tensor_index);
        }

        AxisymmetricConductivityTensors<3> axi_tensors;
        axi_tensors.SetNonConstantConductivities(&non_constant_conductivities);

        axi_tensors.SetFibreOrientationFile("heart/test/data/SimpleOrthotropic3D.fibres");
        TS_ASSERT_THROWS_ANYTHING(axi_tensors.Init());

        axi_tensors.SetFibreOrientationFile("heart/test/data/SimpleAxisymmetric.fibres");
        axi_tensors.Init();

        for (unsigned tensor_index=0; tensor_index<4; tensor_index++)
        {
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](0,0), 100*tensor_index);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](0,1), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](0,2), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](1,0), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](1,1), 10*tensor_index);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](1,2), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](2,0), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](2,1), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](2,2), 10*tensor_index);
        }

    }

};

#endif /*TESTFIBREORIENTATIONTENSORS_HPP_*/
