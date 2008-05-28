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
            auto_ptr<chaste_parameters_type> params (ChasteParameters("heart/test/data/ChasteParametersSlab.xml"));
            
            simulation_type simulation_params = params->Simulation();
            
            TS_ASSERT_EQUALS(simulation_params.SimulationDuration(), 10.0);

            assert(simulation_params.Mesh().LoadMesh() == NULL);
            assert(simulation_params.Mesh().Slab() != NULL);

            TS_ASSERT_EQUALS(simulation_params.Mesh().Slab()->SlabX(), 4.0);
            TS_ASSERT_EQUALS(simulation_params.Mesh().Slab()->SlabY(), 0.1);
            TS_ASSERT_EQUALS(simulation_params.Mesh().Slab()->SlabZ(), 2.0);            
            TS_ASSERT_EQUALS(simulation_params.Mesh().Slab()->InterNodeSpace(), 0.1);

			physiological_type physiological_params = params->Physiological();

            TS_ASSERT_EQUALS(physiological_params.IntracellularConductivities().longi(), 1.75);
            TS_ASSERT_EQUALS(physiological_params.IntracellularConductivities().trans(), 1.75);
            TS_ASSERT_EQUALS(physiological_params.IntracellularConductivities().normal(), 1.75);

            TS_ASSERT_EQUALS(physiological_params.ExtracellularConductivities().longi(), 7.0);
            TS_ASSERT_EQUALS(physiological_params.ExtracellularConductivities().trans(), 7.0);
            TS_ASSERT_EQUALS(physiological_params.ExtracellularConductivities().normal(), 7.0);

            TS_ASSERT_EQUALS(simulation_params.OutputDirectory(), "ChasteResults");
            TS_ASSERT_EQUALS(simulation_params.MeshOutputDirectory(), "Slab");
            TS_ASSERT_EQUALS(simulation_params.Domain(), domain_type::Mono);
            TS_ASSERT_EQUALS(simulation_params.IonicModel(), ionic_model_type::FaberRudy2000Version3);
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
            auto_ptr<chaste_parameters_type> params (ChasteParameters("heart/test/data/ChasteParametersLoadMesh.xml"));
            
            simulation_type simulation_params = params->Simulation();
            
            TS_ASSERT_EQUALS(simulation_params.SimulationDuration(), 10.0);

            assert(simulation_params.Mesh().LoadMesh() != NULL);
            assert(simulation_params.Mesh().Slab() == NULL);

            TS_ASSERT_EQUALS(simulation_params.Mesh().LoadMesh()->name(), "foo");
            TS_ASSERT_EQUALS(simulation_params.Mesh().LoadMesh()->media(), "Orthotropic"); // Testing for the default value
            
            physiological_type physiological_params = params->Physiological();
            
            TS_ASSERT_EQUALS(physiological_params.IntracellularConductivities().longi(), 1.75);
            TS_ASSERT_EQUALS(physiological_params.IntracellularConductivities().trans(), 1.75);
            TS_ASSERT_EQUALS(physiological_params.IntracellularConductivities().normal(), 1.75);

            TS_ASSERT_EQUALS(physiological_params.ExtracellularConductivities().longi(), 7.0);
            TS_ASSERT_EQUALS(physiological_params.ExtracellularConductivities().trans(), 7.0);
            TS_ASSERT_EQUALS(physiological_params.ExtracellularConductivities().normal(), 7.0);

            TS_ASSERT_EQUALS(simulation_params.OutputDirectory(), "ChasteResults");
            TS_ASSERT_EQUALS(simulation_params.MeshOutputDirectory(), "Slab");
            TS_ASSERT_EQUALS(simulation_params.Domain(), domain_type::Mono);
            TS_ASSERT_EQUALS(simulation_params.IonicModel(), ionic_model_type::FaberRudy2000Version3);
        }
        catch (const xml_schema::exception& e)
        {
            std::cerr << e << std::endl;
            TS_FAIL("Schema exception");
        }
    }

    void TestUpdate()
    {
        try
        {
            auto_ptr<chaste_parameters_type> params (ChasteParameters("heart/test/data/ChasteParametersLoadMesh.xml"));
            physiological_type physiological_params = params->Physiological();

            physiological_params.ExtracellularConductivities().longi() = 9.0;
            TS_ASSERT_EQUALS(physiological_params.ExtracellularConductivities().longi(), 9.0);
        }
        catch (const xml_schema::exception& e)
        {
            std::cerr << e << std::endl;
            TS_FAIL("Schema exception");
        }
    }

};
#endif /*TESTSPIRALPARAMETERSREADER_HPP_*/
