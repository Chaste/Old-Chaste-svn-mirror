/*
Copyright (C) University of Oxford, 2008

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with Chaste.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef TESTSPIRALPARAMETERSREADER_HPP_
#define TESTSPIRALPARAMETERSREADER_HPP_

#include <cxxtest/TestSuite.h>
#include <memory>
#include "ChasteParameters.hpp"


using std::auto_ptr;

class TestChasteParametersReader : public CxxTest::TestSuite
{
public:
    void TestReadWithSlab()
    {
        try
        {
            auto_ptr<chaste_parameters_type> p (ChasteParameters("heart/test/data/ChasteParametersSlab.xml"));
            
            TS_ASSERT_EQUALS(p->SimulationDuration(), 10.0);

            assert(p->Mesh().LoadMesh() == NULL);
            assert(p->Mesh().Slab() != NULL);

            TS_ASSERT_EQUALS(p->Mesh().Slab()->SlabX(), 4.0);
            TS_ASSERT_EQUALS(p->Mesh().Slab()->SlabY(), 0.1);
            TS_ASSERT_EQUALS(p->Mesh().Slab()->SlabZ(), 2.0);            
            TS_ASSERT_EQUALS(p->Mesh().Slab()->InterNodeSpace(), 0.1);

            TS_ASSERT_EQUALS(p->IntracellularConductivities().longi(), 1.75);
            TS_ASSERT_EQUALS(p->IntracellularConductivities().trans(), 1.75);
            TS_ASSERT_EQUALS(p->IntracellularConductivities().normal(), 1.75);

            TS_ASSERT_EQUALS(p->ExtracellularConductivities().longi(), 7.0);
            TS_ASSERT_EQUALS(p->ExtracellularConductivities().trans(), 7.0);
            TS_ASSERT_EQUALS(p->ExtracellularConductivities().normal(), 7.0);

            TS_ASSERT_EQUALS(p->OutputDirectory(), "ChasteResults");
            TS_ASSERT_EQUALS(p->MeshOutputDirectory(), "Slab");
            TS_ASSERT_EQUALS(p->Domain(), domain_type::Mono);
            TS_ASSERT_EQUALS(p->IonicModel(), ionic_model_type::FaberRudy2000Version3);
        }
        catch (const xml_schema::exception& e)
        {
            std::cerr << e << std::endl;
            TS_FAIL("Schema exception");
        }
    }

    void TestReadWithLoadMesh()
    {
        try
        {
            auto_ptr<chaste_parameters_type> p (ChasteParameters("heart/test/data/ChasteParametersLoadMesh.xml"));
            
            TS_ASSERT_EQUALS(p->SimulationDuration(), 10.0);

            assert(p->Mesh().LoadMesh() != NULL);
            assert(p->Mesh().Slab() == NULL);

            TS_ASSERT_EQUALS(p->Mesh().LoadMesh()->name(), "foo");
            TS_ASSERT_EQUALS(p->Mesh().LoadMesh()->media(), "Orthotropic"); // Testing for the default value
            
            TS_ASSERT_EQUALS(p->IntracellularConductivities().longi(), 1.75);
            TS_ASSERT_EQUALS(p->IntracellularConductivities().trans(), 1.75);
            TS_ASSERT_EQUALS(p->IntracellularConductivities().normal(), 1.75);

            TS_ASSERT_EQUALS(p->ExtracellularConductivities().longi(), 7.0);
            TS_ASSERT_EQUALS(p->ExtracellularConductivities().trans(), 7.0);
            TS_ASSERT_EQUALS(p->ExtracellularConductivities().normal(), 7.0);

            TS_ASSERT_EQUALS(p->OutputDirectory(), "ChasteResults");
            TS_ASSERT_EQUALS(p->MeshOutputDirectory(), "Slab");
            TS_ASSERT_EQUALS(p->Domain(), domain_type::Mono);
            TS_ASSERT_EQUALS(p->IonicModel(), ionic_model_type::FaberRudy2000Version3);
        }
        catch (const xml_schema::exception& e)
        {
            std::cerr << e << std::endl;
            TS_FAIL("Schema exception");
        }
    }

};
#endif /*TESTSPIRALPARAMETERSREADER_HPP_*/
